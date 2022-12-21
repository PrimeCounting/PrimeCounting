# Prime Counting

This is a repository for the [paper](https://arxiv.org/abs/2212.09857) Computing $\pi(N)$: An elementary approach in $\tilde{O}(\sqrt{N})$ time by Dean Hirsch, Ido Kessler and Uri Mendlovic.

## Introduction

We present an efficient and elementary algorithm for computing the number of primes up to $N$ in $\tilde{O}(\sqrt N)$ time, improving upon the existing combinatorial methods that require $\tilde{O}(N ^ {2/3})$ time. Our method has a similar complexity to the analytical approach to prime counting, while avoiding complex analysis and the use of arbitrary precision complex numbers. 

## Content

The repository contain two implementations of the algorithm described in the paper. One implementation is in python, which is provided for ease of reading, and one in c++ which provides more optimizations.

It is recommended reading the section [`Basic Algorithm`](https://arxiv.org/abs/2212.09857) in the paper before reading the code.

## Build instructions

Requires cmake and a C++ compiler.
```bash
mkdir build
cd build
cmake ../cpp -DCMAKE_CXX_COMPILER=g++-10
make -j countprimes

# To run tests:
cmake ../cpp -DCMAKE_CXX_COMPILER=g++-10 -DCOMPILE_TESTS=ON
make -j tests
```

## Usage
```bash
# Count primes up to 1e10 
./countprimes 1e10
# 10000000000
# Ntt of primes: 100% [0 (s)] 
# Mobius (1025, 100000): 100% [0 (s)] 
# INTT finalize mobius: 100% [0 (s)] 
# SmallPrimeNaiveConvolution: 100% [1 (s)] 
# Error correction: 100% [0 (s)] 
# Num primes up to 10000000000:
#         455052511

# Count primes up to 1e14, use x5 less memory (traded for runtime)
./countprimes 1e14 5
# 100000000000000
# Ntt of primes: 100% [11 (s)] 
# Mobius (1025, 10000000): 100% [8 (s)] 
# INTT finalize mobius: 100% [11 (s)] 
# SmallPrimeNaiveConvolution: 100% [23 (s)] 
# Error correction: 100% [91 (s)] 
# Num primes up to 100000000000000:
#     3204941750802
```

## Benchmark
| Up To       | Time        |
| ----------- | ----------- |
| 1e7         | 0.06 sec    |
| 1e8         | 0.14 sec    |
| 1e9         | 0.35 sec    |
| 1e10        | 1.23 sec    |
| 1e11        | 4.01 sec    |
| 1e12        | 14.6 sec    |
| 1e13        | 46.1 sec    |
| 1e14        | 127 sec     |

Times were taken using a single core on `Intel(R) Core(TM) i3-1005G1 CPU @ 1.20GHz`, with memory tradeoff set to `5`.
