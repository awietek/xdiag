<!--
SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>

SPDX-License-Identifier: Apache-2.0
-->

![cpp](https://img.shields.io/badge/C++-17-blue.svg?style=flat&logo=c%2B%2B)
[![Linux CI](https://github.com/awietek/xdiag/actions/workflows/linux.yml/badge.svg?style=for-the-badge)](https://github.com/awietek/xdiag/actions/workflows/linux.yml)
[![Mac OSX CI](https://github.com/awietek/xdiag/actions/workflows/osx.yml/badge.svg?style=for-the-badge)](https://github.com/awietek/xdiag/actions/workflows/osx.yml)
[![Intel MPI CI](https://github.com/awietek/xdiag/actions/workflows/intelmpi.yml/badge.svg?style=for-the-badge)](https://github.com/awietek/xdiag/actions/workflows/intelmpi.yml)
[![Examples CI](https://github.com/awietek/xdiag/actions/workflows/examples.yml/badge.svg?style=for-the-badge)](https://github.com/awietek/xdiag/actions/workflows/examples.yml)

# XDiag
## High-performance Yxact Diagonalization Routines and Algorithms

A C++ library to perform efficient Exact Diagonalizations of quantum many body systems. 


| **Documentation**                                                                          | **Citation**                                                                                           | **License**                                                        |
|:------------------------------------------------------------------------------------------:|:------------------------------------------------------------------------------------------------------:|--------------------------------------------------------------------|
| [![docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://awietek.github.io/xdiag) | [![arXiv](https://img.shields.io/badge/arXiv-2007.14822-b31b1b.svg)](https://arxiv.org/abs/2505.02901) | ![license](https://img.shields.io/badge/license-Apache%202.0-blue) |


### Features:
- Basic algebra of operators in quantum many-body systems
- Iterative linear algebra for computing eigendecompositions and time-evolutions (e.g. Lanczos algorithm)
- Local spin, t-J, or fermionic models
- Full support of generic space group symmetries
- parallelization both with OpenMP and MPI
- modern C++17 impementation simplifying usage
- wrapped in a convenient Julia library: [XDiag.jl](https://github.com/awietek/XDiag.jl)

### Installation:
Clone this repository first. Afterwards, the **xdiag** library can be compiled using the standard CMake instructions
```bash
cmake -S . -B build
cmake --build build
cmake --install build
```

### Example Code:
```cpp
#include <xdiag/all.hpp>

using namespace xdiag;

int main() try {
  
  int nsites = 16;
  int nup = nsites / 2;
  Spinhalf block(nsites, nup);

  // Define the nearest-neighbor Heisenberg model
  OpSum ops;
  for (int i = 0; i < nsites; ++i) {
    ops += "J" * Op("SdotS", {i, (i + 1) % nsites});
  }
  ops["J"] = 1.0;

  set_verbosity(2);                // set verbosity for monitoring progress
  double e0 = eigval0(ops, block); // compute ground state energy
  
  Log("Ground state energy: {:.12f}", e0);
  
} catch (Error e) {
  error_trace(e);
}

```

### Documentation
The full documentation is available at [awietek.github.io/xdiag](https://awietek.github.io/xdiag).

### About
author:   Alexander Wietek
license:   Apache License 2.0
