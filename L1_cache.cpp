#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <string>
#include <queue>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <limits>
#include <getopt.h>

// MESI Protocol states
enum class MESIState {
    INVALID,
    SHARED,
    EXCLUSIVE,
    MODIFIED
};

// Memory operation type
enum class MemOperation {
    READ,
    WRITE
};

// Bus operation type
enum class BusOperation {
    BUSRD,        // Bus Read
    BUSRDX,       // Bus Read Exclusive
    BUSUPGR,      // Bus Upgrade
};

// Response type (for snooping)
enum class SnoopResult {
    NONE,           // No response needed
    SHARED,         // Block is shared
    MODIFIED,       // Block is modified (needs to provide data)
    EXCLUSIVE       // Block is exclusive (optional response)
};

// Cache line structure
struct CacheLine {
    uint32_t tag;
    MESIState state;
    std::vector<uint8_t> data;
    uint64_t lastUsed;

    explicit CacheLine(const int blockSize) : tag(0), state(MESIState::INVALID), data(blockSize, 0), lastUsed(0) {}
};

// Cache set structure
struct CacheSet {
    std::vector<CacheLine> lines;

    CacheSet(const int E, int blockSize) {
        for (int i = 0; i < E; i++) {
            lines.emplace_back(blockSize);
        }
    }
};

// Forward declaration
class CacheSimulator;

// Cache structure
class Cache {
    int S; // Number of sets (2^s)
    int E; // Associativity
    int B; // Block size (2^b)
    int s; // Number of set index bits
    int b; // Number of block offset bits
    int t; // Number of tag bits

    std::vector<CacheSet> sets;

    // Cache statistics
    int coreId;
    uint64_t readCount;
    uint64_t writeCount;
    uint64_t missCount;
    uint64_t hitCount;
    uint64_t evictionCount;
    uint64_t writebackCount;
    uint64_t cyclesCount;
    uint64_t idleCycles;
    uint64_t busInvalidations;
    uint64_t busDataTraffic;
    uint64_t busTransactions;

public:
    Cache(const int coreId, const int s, int E, const int b) : E(E), s(s), b(b), coreId(coreId) {
        S = 1 << s;
        B = 1 << b;
        t = 32 - s - b;

        sets.reserve(S);
        for (int i = 0; i < S; i++) {
            sets.emplace_back(E, B);
        }

        // Initialize statistics
        readCount = 0;
        writeCount = 0;
        missCount = 0;
        hitCount = 0;
        evictionCount = 0;
        writebackCount = 0;
        cyclesCount = 0;
        idleCycles = 0;
        busInvalidations = 0;
        busDataTraffic = 0;
        busTransactions = 0;
    }

    // Parse address to get tag, set index, and block offset
    void parseAddress(const uint32_t address, uint32_t& tag, uint32_t& setIndex, uint32_t& blockOffset) const {
        blockOffset = address & (1 << b) - 1;
        setIndex = address >> b & (1 << s) - 1;
        tag = address >> (s + b);
    }

    // Get block address (clear block offset bits)
    [[nodiscard]] uint32_t getBlockAddress(const uint32_t address) const {
        return address & ~((1 << b) - 1);
    }

    // Check if a cache line is present in the cache
    std::pair<bool, int> isHit(const uint32_t address, const uint64_t cycle) {
        uint32_t tag, setIndex, blockOffset;
        parseAddress(address, tag, setIndex, blockOffset);

        CacheSet& set = sets[setIndex];
        for (int i = 0; i < E; i++) {
            if (set.lines[i].tag == tag &&
                set.lines[i].state != MESIState::INVALID) {
                set.lines[i].lastUsed = cycle;
                return {true, i};
            }
        }
        return {false, -1};
    }

    // Find LRU line for replacement
    [[nodiscard]] int findLRULine(const uint32_t setIndex) const {
        const CacheSet& set = sets[setIndex];
        int lru_line = 0;
        uint64_t lruTime = std::numeric_limits<uint64_t>::max();

        for (int i = 0; i < E; i++) {
            if (set.lines[i].state == MESIState::INVALID) {
                return i; // Empty line found
            }
            if (set.lines[i].lastUsed < lruTime) {
                lruTime = set.lines[i].lastUsed;
                lru_line = i;
            }
        }
        return lru_line;
    }

    // Process a memory operation (read or write)
    std::pair<bool, BusOperation> processMemoryOperation(uint32_t address, MemOperation op,
                                                        const std::vector<Cache*>& otherCaches, uint64_t& cycles);

    // Check if an operation needs the bus
    bool needsBus(const uint32_t address, const MemOperation op, const uint64_t cycle) {
        uint32_t tag, setIndex, blockOffset;
        parseAddress(address, tag, setIndex, blockOffset);

        // Check cache hit/miss
        auto [hit, lineIndex] = isHit(address, cycle);

        if (op == MemOperation::READ) {
            if (hit) {
                return false; // Read hit doesn't need bus
            }
            return true; // Read miss needs bus
        }
        // WRITE operation
        if (hit) {
            const CacheLine& line = sets[setIndex].lines[lineIndex];

            if (line.state == MESIState::MODIFIED || line.state == MESIState::EXCLUSIVE) {
                return false; // Write hit to M or E doesn't need bus
            }
            return true; // Write hit to S needs bus for invalidation
        }
        return true; // Write miss needs bus
    }

    // Handle snooping of bus transactions - returns SnoopResult
    std::pair<SnoopResult, bool> handleBusOperation(const uint32_t address, const BusOperation busOp, int requestorId) {
        uint32_t tag, setIndex, blockOffset;
        parseAddress(address, tag, setIndex, blockOffset);

        auto [hit, lineIndex] = isHit(address, 0); // Using 0 as cycle doesn't matter for snooping

        if (!hit) {
            return {SnoopResult::NONE, false}; // Block not in cache, nothing to do
        }

        CacheLine& line = sets[setIndex].lines[lineIndex];
        bool wasInvalidated = false;

        switch (busOp) {
            case BusOperation::BUSRD: {
                if (line.state == MESIState::MODIFIED) {
                    // Provide data and transition to SHARED
                    line.state = MESIState::SHARED;
                    return {SnoopResult::MODIFIED, false};
                }
                if (line.state == MESIState::EXCLUSIVE) {
                    // Transition to SHARED
                    line.state = MESIState::SHARED;
                    return {SnoopResult::EXCLUSIVE, false};
                }
                if (line.state == MESIState::SHARED) {
                    return {SnoopResult::SHARED, false};
                }
            }
            break;

            case BusOperation::BUSRDX:
                if (line.state != MESIState::INVALID) {
                    SnoopResult result = line.state == MESIState::MODIFIED ?
                                           SnoopResult::MODIFIED : SnoopResult::SHARED;

                    if (line.state == MESIState::MODIFIED) {
                        writebackCount++; // Implicit writeback to memory
                    } else if (line.state == MESIState::EXCLUSIVE || line.state == MESIState::SHARED) {
                    }

                    line.state = MESIState::INVALID;
                    wasInvalidated = true;
                    return {result, wasInvalidated};
                }
                break;

            case BusOperation::BUSUPGR:
                if (line.state == MESIState::SHARED) {
                    // Another cache wants to upgrade SHARED->MODIFIED
                    line.state = MESIState::INVALID;
                    wasInvalidated = true;
                    return {SnoopResult::SHARED, wasInvalidated};
                }
                break;

            default:
                break;
        }

        return {SnoopResult::NONE, false};
    }

    // Register an invalidation
    void registerInvalidation() {
        busInvalidations++;
    }

    // Register data traffic
    void registerDataTraffic(const int size) {
        busDataTraffic += size;
    }

    // Increment bus transactions
    void incrementBusTransactions() {
        busTransactions++;
    }

    // Getters for statistics
    [[nodiscard]] uint64_t getReadCount() const { return readCount; }
    [[nodiscard]] uint64_t getWriteCount() const { return writeCount; }
    [[nodiscard]] uint64_t getMissCount() const { return missCount; }
    [[nodiscard]] uint64_t getHitCount() const { return hitCount; }
    [[nodiscard]] uint64_t getEvictionCount() const { return evictionCount; }
    [[nodiscard]] uint64_t getWritebackCount() const { return writebackCount; }
    [[nodiscard]] uint64_t getCyclesCount() const { return cyclesCount; }
    [[nodiscard]] uint64_t getIdleCycles() const { return idleCycles; }
    [[nodiscard]] uint64_t getBusInvalidations() const { return busInvalidations; }
    [[nodiscard]] uint64_t getBusDataTraffic() const { return busDataTraffic; }
    [[nodiscard]] uint64_t getBusTransactions() const { return busTransactions; }

    void updateCyclesCount(const uint64_t cycles) {
        cyclesCount = cycles;
    }

    void incrementIdleCycles() {
        idleCycles++;
    }

    [[nodiscard]] int getCoreId() const { return coreId; }
    [[nodiscard]] int getBlockSize() const { return B; }
};

// Memory trace structure
struct MemoryTrace {
    MemOperation op;
    uint32_t address;
};

// System simulator
class CacheSimulator {
    std::vector<Cache> caches;
    std::vector<std::vector<MemoryTrace>> traces;
    int s; // Number of set index bits
    int E; // Associativity
    int b; // Number of block offset bits

    // Bus state variables
    int busOwner;             // Core ID of the current bus owner (-1 if bus is free)
    uint64_t busFreeAtCycle;  // Cycle when the bus will be free again

public:
    CacheSimulator(int s, int E, int b) : s(s), E(E), b(b),
                                         busOwner(-1), busFreeAtCycle(0) {
        // Create 4 caches (quad-core)
        for (int i = 0; i < 4; i++) {
            caches.emplace_back(i, s, E, b);
        }

        // Initialize traces for 4 cores
        traces.resize(4);
    }

    // Load traces from files
    bool loadTraces(const std::string& baseFileName) {
        for (int i = 0; i < 4; i++) {
            std::string fileName = baseFileName + "_proc" + std::to_string(i) + ".trace";
            std::ifstream file(fileName);

            if (!file.is_open()) {
                std::cerr << "Error: Could not open file " << fileName << std::endl;
                return false;
            }

            std::string line;
            while (std::getline(file, line)) {
                if (line.empty()) continue;

                char op;
                std::string addrStr;
                std::istringstream iss(line);

                iss >> op >> addrStr;

                if (op != 'R' && op != 'W') {
                    std::cerr << "Error: Invalid operation " << op << " in file " << fileName << std::endl;
                    continue;
                }

                MemoryTrace trace{};
                trace.op = op == 'R' ? MemOperation::READ : MemOperation::WRITE;

                // Parse address from hex string
                trace.address = 0;
                if (addrStr.substr(0, 2) == "0x") {
                    addrStr = addrStr.substr(2);
                }

                try {
                    trace.address = std::stoul(addrStr, nullptr, 16);
                } catch (const std::exception& _) {
                    std::cerr << "Error parsing address: " << addrStr << " in file " << fileName << std::endl;
                    continue;
                }

                traces[i].push_back(trace);
            }

            file.close();
        }

        return true;
    }

    // Run simulation
    void runSimulation() {
        std::vector<size_t> positions(4, 0);      // Current position in each trace
        std::vector<uint64_t> cycles(4, 0);       // Current cycle count for each core
        std::vector waiting(4, false);      // Is a core waiting for the bus?
        bool done = false;
        uint64_t globalCycle = 0;

        // Create pointers to other caches for each core
        std::vector<std::vector<Cache*>> otherCaches(4);
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                if (i != j) {
                    otherCaches[i].push_back(&caches[j]);
                }
            }
        }

        while (!done) {
            globalCycle++;

            // Check if the bus is free at this cycle
            if (busOwner != -1 && globalCycle >= busFreeAtCycle) {
                busOwner = -1; // Bus is now free
            }

            // Find cores that are ready to execute an instruction
            std::vector<int> readyCores;
            for (int i = 0; i < 4; i++) {
                if (positions[i] < traces[i].size() && cycles[i] <= globalCycle) {
                    readyCores.push_back(i);
                }
            }

            // Sort by core ID for priority
            std::sort(readyCores.begin(), readyCores.end());

            // Check which cores need the bus
            for (const int coreId : readyCores) {
                if (positions[coreId] >= traces[coreId].size()) {
                    continue; // Skip if already finished
                }

                const auto&[op, address] = traces[coreId][positions[coreId]];
                const bool needsBus = caches[coreId].needsBus(address, op, cycles[coreId]);

                if (!needsBus) {
                    // Process instruction that doesn't need the bus

                    caches[coreId].processMemoryOperation(
                        address, op, otherCaches[coreId], cycles[coreId]);

                    // Advance to next instruction
                    positions[coreId]++;

                    // Update cache cycle count
                    caches[coreId].updateCyclesCount(cycles[coreId]);

                } else if (busOwner == -1) {
                    // This core gets the bus
                    busOwner = coreId;

                    auto [usedBus, busOp] = caches[coreId].processMemoryOperation(
                        address, op, otherCaches[coreId], cycles[coreId]);

                    // Handle snooping on other caches for block address (not the specific word)
                    const uint32_t blockAddress = caches[coreId].getBlockAddress(address);

                    // Track invalidations performed by this operation
                    int invalidationsDone = 0;

                    for (Cache* otherCache : otherCaches[coreId]) {
                        auto [result, wasInvalidated] = otherCache->handleBusOperation(blockAddress, busOp, coreId);

                        if (wasInvalidated) {
                            invalidationsDone++;
                        }
                    }

                    // Register invalidations in the requesting core
                    if (invalidationsDone > 0) {
                        caches[coreId].registerInvalidation();
                    }

                    // Bus is busy until this operation completes
                    busFreeAtCycle = cycles[coreId];

                    // Advance to next instruction
                    positions[coreId]++;

                    // Update cache cycle count
                    caches[coreId].updateCyclesCount(cycles[coreId]);

                } else {
                    // This core needs the bus but can't get it - it's idle
                    caches[coreId].incrementIdleCycles();
                    cycles[coreId]++;
                }
            }

            // Check if all cores have finished
            done = true;
            for (int i = 0; i < 4; i++) {
                if (positions[i] < traces[i].size()) {
                    done = false;
                    break;
                }
            }
        }
    }

    // Print simulation results
    void printResults(std::ostream& out) const {
        // Calculate cache size in KB per core
        const double cacheSize = (1 << s) * E * (1 << b) / 1024.0;

        // Calculate total bus transactions and traffic
        uint64_t totalBusTransactions = 0;
        uint64_t totalBusTraffic = 0;
        for (int i = 0; i < 4; i++) {
            totalBusTransactions += caches[i].getBusTransactions();
            totalBusTraffic += caches[i].getBusDataTraffic();
        }

        // Print simulation parameters
        out << "Simulation Parameters:" << std::endl;
        const std::string tracePrefix = "Unknown"; // This would need to be passed in or stored
        out << "Trace Prefix: " << tracePrefix << std::endl;
        out << "Set Index Bits: " << s << std::endl;
        out << "Associativity: " << E << std::endl;
        out << "Block Bits: " << b << std::endl;
        out << "Block Size (Bytes): " << (1 << b) << std::endl;
        out << "Number of Sets: " << (1 << s) << std::endl;
        out << "Cache Size (KB per core): " << std::fixed << std::setprecision(2) << cacheSize << std::endl;
        out << "MESI Protocol: Enabled" << std::endl;
        out << "Write Policy: Write-back, Write-allocate" << std::endl;
        out << "Replacement Policy: LRU" << std::endl;
        out << "Bus: Central snooping bus" << std::endl;
        out << std::endl;

        // Print per-core statistics
        for (int i = 0; i < 4; i++) {
            out << "Core " << i << " Statistics:" << std::endl;
            const uint64_t totalInstructions = caches[i].getReadCount() + caches[i].getWriteCount();
            const uint64_t executionCycles = caches[i].getCyclesCount() - caches[i].getIdleCycles();
            const double missRate = static_cast<double>(caches[i].getMissCount()) / static_cast<double>(totalInstructions) * 100.0;

            out << "Total Instructions: " << totalInstructions << std::endl;
            out << "Total Reads: " << caches[i].getReadCount() << std::endl;
            out << "Total Writes: " << caches[i].getWriteCount() << std::endl;
            out << "Total Execution Cycles: " << executionCycles << std::endl;
            out << "Idle Cycles: " << caches[i].getIdleCycles() << std::endl;
            out << "Cache Misses: " << caches[i].getMissCount() << std::endl;
            out << "Cache Miss Rate: " << std::fixed << std::setprecision(2) << missRate << "%" << std::endl;
            out << "Cache Evictions: " << caches[i].getEvictionCount() << std::endl;
            out << "Writebacks: " << caches[i].getWritebackCount() << std::endl;
            out << "Bus Invalidations: " << caches[i].getBusInvalidations() << std::endl;
            out << "Data Traffic (Bytes): " << caches[i].getBusDataTraffic() << std::endl;
            out << std::endl;
        }

        // Print overall bus summary
        out << "Overall Bus Summary:" << std::endl;
        out << "Total Bus Transactions: " << totalBusTransactions << std::endl;
        out << "Total Bus Traffic (Bytes): " << totalBusTraffic << std::endl;
    }

    // Get maximum execution time (for experiments)
    [[nodiscard]] uint64_t getMaxExecutionTime() const {
        uint64_t maxCycles = 0;
        for (int i = 0; i < 4; i++) {
            maxCycles = std::max(maxCycles, caches[i].getCyclesCount());
        }
        return maxCycles;
    }
};

// Implementation of processMemoryOperation
std::pair<bool, BusOperation> Cache::processMemoryOperation(const uint32_t address, const MemOperation op,
                                                          const std::vector<Cache*>& otherCaches, uint64_t& cycles) {
    uint32_t tag, setIndex, blockOffset;
    parseAddress(address, tag, setIndex, blockOffset);

    // Check cache hit/miss
    auto [hit, lineIndex] = isHit(address, cycles);

    if (op == MemOperation::READ) {
        readCount++;
        if (hit) {
            hitCount++;
            cycles += 1; // Cache query hit
            return {false, BusOperation::BUSRD}; // No bus operation needed
        }
        missCount++;
        // Need to get block from other caches or memory
        cycles += 1; // Cache lookup (checking if block exists in any other cache line) + BusRD (Assumed to occur in same clock cycle)
        incrementBusTransactions();
        // Issue BusRd
        auto busOp = BusOperation::BUSRD;

        // Check if other caches have the block
        bool foundInOtherCache = false;
        bool foundModified = false;
        const int wordsInBlock = B / 4; // 4 bytes per word

        for (Cache* otherCache : otherCaches) {
            auto [otherHit, otherLineIndex] = otherCache->isHit(getBlockAddress(address), cycles);
            if (otherHit) {
                foundInOtherCache = true;

                uint32_t otherTag, otherSetIndex, otherBlockOffset;
                otherCache->parseAddress(address, otherTag, otherSetIndex, otherBlockOffset);

                const MESIState& otherState = otherCache->sets[otherSetIndex].lines[otherLineIndex].state;

                if (otherState == MESIState::MODIFIED) {
                    otherCache->registerDataTraffic(B); // Count data written to RAM
                    foundModified = true;
                }

                break;
            }
        }

        // If found in another cache (modified or not)
        if (foundInOtherCache) {
            registerDataTraffic(B); // Count receiving data from other cache
            // Modified data is provided by other cache, 2N cycles for transfer
            cycles += 2 * wordsInBlock;
            if (foundModified) {
                // Add 100 cycles for memory write back
                writebackCount++;
                cycles += 100;
            }
        } else {
            // Get data from memory
            registerDataTraffic(B); // Count receiving data from memory
            cycles += 100; // Memory access
        }

        // Add the block to cache
        const int lruLine = findLRULine(setIndex);
        CacheLine& line = sets[setIndex].lines[lruLine];

        if (line.state != MESIState::INVALID) {
            // Evict the line
            evictionCount++;

            if (line.state == MESIState::MODIFIED) {
                // Writeback to memory
                writebackCount++;
                registerDataTraffic(B); // Count writing data to memory
                cycles += 100; // Writeback to memory
            }
        }
        line.tag = tag;
        line.lastUsed = cycles;
        line.state = foundInOtherCache ? MESIState::SHARED : MESIState::EXCLUSIVE;

        return {true, busOp};
    }
    // WRITE operation
    writeCount++;

    if (hit) {
        CacheLine& line = sets[setIndex].lines[lineIndex];
        cycles += 1; // Cache hit: 1 cycle
        if (line.state == MESIState::MODIFIED || line.state == MESIState::EXCLUSIVE) {
            // Already modified, just update
            line.state = MESIState::MODIFIED;
            hitCount++;
            return {false, BusOperation::BUSRDX}; // No bus operation needed
        }
        // SHARED state
        // Need to invalidate copies in other caches
        line.state = MESIState::MODIFIED;
        hitCount++;

        // Issue BusUPGR

        auto busOp = BusOperation::BUSUPGR;
        incrementBusTransactions();

        return {true, busOp};
    }
    // WRITE MISS
    missCount++;
    cycles += 1; // 1 cycle for cache lookup + BusRDX

    // Issue BusRdX
    auto busOp = BusOperation::BUSRDX;
    incrementBusTransactions();
    // Check if other caches have the block (for invalidation)
    bool foundInOtherCache = false;
    bool foundModified = false;

    for (Cache* otherCache : otherCaches) {
        auto [otherHit, otherLineIndex] = otherCache->isHit(getBlockAddress(address), cycles);
        if (otherHit) {
            foundInOtherCache = true;

            uint32_t otherTag, otherSetIndex, otherBlockOffset;
            otherCache->parseAddress(address, otherTag, otherSetIndex, otherBlockOffset);

            const MESIState& otherState = otherCache->sets[otherSetIndex].lines[otherLineIndex].state;

            if (otherState == MESIState::MODIFIED) {
                foundModified = true;
                otherCache->registerDataTraffic(B); // count write from other cache into memory
            }

            break;
        }
    }

    // If found in another cache
    if (foundInOtherCache) {
        registerDataTraffic(B); // Count receiving data from memory
        if (foundModified) {
            cycles += 200; // Memory access + Memory write
        } else {
            cycles += 100; // Memory access
        }
    } else {
        // Get data from memory
        registerDataTraffic(B); // Count receiving data from memory
        cycles += 100; // Memory access
    }

    // Add the block to cache
    const int lruLine = findLRULine(setIndex);
    CacheLine& line = sets[setIndex].lines[lruLine];

    if (line.state != MESIState::INVALID) {
        // Evict the line
        evictionCount++;

        if (line.state == MESIState::MODIFIED) {
            // Writeback to memory
            writebackCount++;
            registerDataTraffic(B); // Count writing data to memory
            cycles += 100; // Writeback to memory
        }
    }

    line.tag = tag;
    line.lastUsed = cycles;
    line.state = MESIState::MODIFIED; // Write puts it directly in modified state

    return {true, busOp};
}

void printHelp() {
    std::cout << "Usage: L1simulate [OPTIONS]" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -t <tracefile>: name of parallel application (e.g. app1) whose 4 traces are to be used in simulation" << std::endl;
    std::cout << "  -s <s>: number of set index bits (number of sets in the cache = S = 2^s)" << std::endl;
    std::cout << "  -E <E>: associativity (number of cache lines per set)" << std::endl;
    std::cout << "  -b <b>: number of block bits (block size = B = 2^b)" << std::endl;
    std::cout << "  -o <outfilename>: logs output in file for plotting etc." << std::endl;
    std::cout << "  -h: prints this help" << std::endl;
}

int main(int argc, char* argv[]) {
    std::string traceFile;
    std::string outputFile;
    int s = 6; // Default: 64 sets (2^6)
    int E = 2; // Default: 2-way associative
    int b = 5; // Default: 32-byte block size (2^5)

    // Parse command line arguments
    int opt;
    while ((opt = getopt(argc, argv, "t:s:E:b:o:h")) != -1) {
        switch (opt) {
            case 't':
                traceFile = optarg;
                break;
            case 's':
                s = std::stoi(optarg);
                break;
            case 'E':
                E = std::stoi(optarg);
                break;
            case 'b':
                b = std::stoi(optarg);
                break;
            case 'o':
                outputFile = optarg;
                break;
            case 'h':
                printHelp();
                return 0;
            default:
                std::cerr << "Invalid option" << std::endl;
                printHelp();
                return 1;
        }
    }

    if (traceFile.empty()) {
        std::cerr << "Error: Trace file not specified" << std::endl;
        printHelp();
        return 1;
    }

    // Create cache simulator
    CacheSimulator simulator(s, E, b);

    // Load traces
    if (!simulator.loadTraces(traceFile)) {
        std::cerr << "Error loading traces" << std::endl;
        return 1;
    }

    // Run simulation
    simulator.runSimulation();

    // Print results
    if (!outputFile.empty()) {
        std::ofstream outFile(outputFile);
        if (outFile.is_open()) {
            simulator.printResults(outFile);
            outFile.close();
        } else {
            std::cerr << "Error: Could not open output file " << outputFile << std::endl;
            simulator.printResults(std::cout);
        }
    } else {
        simulator.printResults(std::cout);
    }

    return 0;
}