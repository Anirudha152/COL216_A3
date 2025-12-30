# COL216 Assignment 3: Multicore Cache Simulator with MESI Coherence Protocol

**Authors:** Anirudha Saraf, Ayush Kumar Singh

This project implements a multicore cache simulator with MESI coherence protocol support for analyzing cache performance in parallel systems.

## Overview

The L1 Cache Simulator models a quad-core system with private L1 caches connected via a shared bus. It implements:
- **MESI coherence protocol** for maintaining cache coherence
- **LRU replacement policy** for cache line eviction
- **Detailed statistics tracking** for performance analysis
- **Parameterized cache geometry** (sets, associativity, block size)
- **Memory traces** for realistic workload simulation

## Features

- Accurately simulates the behavior of a 4-core system with private L1 caches
- Implements the MESI (Modified, Exclusive, Shared, Invalid) coherence protocol
- Models bus contention and coherence traffic
- Tracks comprehensive statistics including:
    - Cache hit/miss rates
    - Execution and idle cycles
    - Bus transactions and data traffic
    - Cache evictions and writebacks
- Supports configurable cache parameters (size, associativity, block size)
- Processes memory traces to simulate realistic application behavior

## Requirements

- C++17 compatible compiler (g++ recommended)
- Standard C++ libraries

## Building

To build the simulator, simply run:

```bash
make
```

This will compile the source code and create the `L1simulate` executable.

To clean build artifacts:

```bash
make clean
```

## Usage

```
Usage: L1simulate [OPTIONS]
Options:
  -t <tracefile>: name of parallel application (e.g. app1) whose 4 traces are to be used in simulation
  -s <s>: number of set index bits (number of sets in the cache = S = 2^s)
  -E <E>: associativity (number of cache lines per set)
  -b <b>: number of block bits (block size = B = 2^b)
  -o <outfilename>: logs output in file for plotting etc.
  -h: prints this help
```

### Trace File Format

The simulator reads trace files named `<tracefile>_proc<N>.trace` where N ranges from 0 to 3 (one file per core). Each line of a trace file should have the format:

```
<operation> <address>
```

Where:
- `<operation>` is either `R` (for read) or `W` (for write)
- `<address>` is the memory address in hexadecimal (with or without 0x prefix)

Example:
```
R 0x1a2b3c4d
W 0xdeadbeef
```

## Example

To run the simulator with a trace file named "matmul" with 64 sets (s=6), 2-way associativity (E=2), and 32-byte blocks (b=5):

```bash
./L1simulate -t matmul -s 6 -E 2 -b 5 -o results.txt
```

This will simulate the execution of the matmul application on four cores and write the results to results.txt.

## Output Format

The simulator outputs detailed statistics about the simulation, including:

- Simulation parameters (cache configuration, protocols, etc.)
- Per-core statistics (instructions, cycles, miss rates, etc.)
- Overall bus summary (transactions, data traffic)

## Design and Implementation Details

### Cache Structure
- Each cache is configured with parameters s (set index bits), E (associativity), and b (block offset bits)
- Number of sets S = 2^s
- Block size B = 2^b bytes
- Total cache size = S × E × B bytes

### MESI Protocol
The simulator implements the full MESI protocol with four states:
- **M (Modified)**: Cache line is modified and not in other caches
- **E (Exclusive)**: Cache line is unmodified and not in other caches
- **S (Shared)**: Cache line is unmodified and may exist in other caches
- **I (Invalid)**: Cache line is invalid or does not contain valid data

### Bus Operations
- **BusRD**: Bus read for read misses
- **BusRDX**: Bus read exclusive for write misses
- **BusUPGR**: Bus upgrade for write hits to shared blocks

### Timing Assumptions
- Read/Write to memory: 100 cycles
- Cache-to-cache transfer: 2 cycles per word (4 bytes)
- Read Hit: 1 cycle (No bus transactions)
- Write Hit to M/E: 1 cycle (No bus transactions)
- Write Hit to S: 1 cycle + bus operations
- Invalidation: 0 cycles (asynchronous process)
- Lookup + Bus Signal: 1 cycle

## Disclaimer
This code was developed as part of an academic assignment for the COL216 course at IIT Delhi. It is intended for educational purposes and was completed in May 2025, it is no longer being actively updated or maintained, please reach out to me over email or linkedin for any queries. Please cite appropriately if used in research or projects.
