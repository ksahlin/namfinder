# namfinder: A fast NAM finder based on strobemers

Namfinder is mapping tool used to find non-overlapping approximate matches (NAMs). NAMs are constructed from overlapping strobemer seeds. 

Namfinder is a more optimized version of the previous proof-of-concept tool StrobeMap that was implemented for the strobemers paper.

Namfinder has borrowed the whole indexing construction codebase from [strobealign](https://github.com/ksahlin/strobealign) (a short-read mapper), but is used only for finding NAM seeds. 
Credits to @marcelm, @luispedro and @psj1997 for the optimized indexing implementation.


## Features

- Multithreading support
- Fast indexing (2-5 minutes for a human-sized reference genome)
- Output in MUMmer MEM tsv format


## Table of contents

1. [Installation](#installation)
2. [Usage](#usage)
3. [Command-line options](#command-line-options)
4. [Index file](#index-files)
5. [Changelog](#changelog)
6. [Contributing](#contributing)
7. [Performance](#v07-performance)
8. [Credits](#credits)
9. [Version info](#version-info)
10. [License](#licence)

## Installation

You need to have CMake, a recent `g++` (tested with version 8) and [zlib](https://zlib.net/) installed.
Then do the following:
```
git clone https://github.com/ksahlin/namfinder
cd namfinder
cmake -B build -DCMAKE_C_FLAGS="-march=native" -DCMAKE_CXX_FLAGS="-march=native"
make -j -C build
```
The resulting binary is `build/namfinder`.

The binary is tailored to the CPU the compiler runs on.
If it needs to run on other machines, use this `cmake` command instead for compatibility with most x86-64 CPUs in use today:
```
cmake -B build -DCMAKE_C_FLAGS="-msse4.2" -DCMAKE_CXX_FLAGS="-msse4.2"
```


## Usage

Parameter `-k` is the strobe size, `-s` is sub-k-mer size (used for thinning in syncmers). Set `-s` to the same value as `k`for no thinning.
Parameters `-l` and `-u` are window min and window mac for sampling the downstream strobe. only strobemers of order 2 can currently be used.


```
namfinder -k 10 -s 10 -l 11 -u 35 -C 500 -o nams.tsv ref.fa reads.f[a/q]
```



CREDITS
----------------


- Some of the ideas for the index and NAM 
construction in namfinder was borrowed from: 
_Sahlin, K. Strobealign: flexible seed size enables ultra-fast and accurate read alignment. 
Genome Biol 23, 260 (2022)._ 
[https://doi.org/10.1186/s13059-022-02831-7](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02831-7)
- Big improvements were designed by @marcelm and @luispedro, and inplemented by @marcelm and @psj1997 (forthcoming paper). 

LICENCE
----------------

MIT license, see [LICENSE](https://github.com/ksahlin/strobealign/blob/main/LICENSE).

