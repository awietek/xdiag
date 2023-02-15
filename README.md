![license](https://img.shields.io/badge/license-Apache%202.0-blue)
![cpp](https://img.shields.io/badge/C++-17-blue.svg?style=flat&logo=c%2B%2B)
[![docs](https://img.shields.io/badge/Documentation-here-red.svg)](https://awietek.github.io/hydradoc)
[![Linux CI](https://github.com/awietek/hydra/actions/workflows/linux.yml/badge.svg?style=for-the-badge)](https://github.com/awietek/hydra/actions/workflows/linux.yml)
[![Mac OSX CI](https://github.com/awietek/hydra/actions/workflows/osx.yml/badge.svg?style=for-the-badge)](https://github.com/awietek/hydra/actions/workflows/osx.yml)
[![DOI](https://zenodo.org/badge/169422780.svg)](https://zenodo.org/badge/latestdoi/169422780)


# Hydra
## High-performance Yxact Diagonalization Routines and Algorithms

A C++ library to perform efficient Exact Diagonalizations of quantum many body systems. 

### Features:
- Basic algebra of operators in quantum many-body systems
- Iterative linear algebra for computing eigendecompositions and time-evolutions (e.g. Lanczos algorithm)
- Local spin, t-J, or fermionic models
- Full support of generic space group symmetries
- parallelization both with OpenMP and MPI
- modern C++17 impementation simplifying usage

### Installation:
Clone this repository first. Afterwards, the **hydra** library can be compiled using the standard CMake instructions
```bash
mkdir build
cd build
cmake ..
make
```

### Example Code:
```cpp
#include <hydra/all.h>

int main() {
  using namespace hydra;

  int n_sites = 16;
  int nup = n_sites / 2;

  // Define the Hilbert space block
  auto block = Spinhalf(n_sites, nup);

  // Define the nearest-neighbor Heisenberg Hamiltonian
  BondList bonds;
  for (int i = 0; i < n_sites; ++i) {
    bonds << Bond("HB", "J", {i, (i + 1) % n_sites});
  }

  // Set the coupling constant "J" to one
  bonds["J"] = 1.0;

  // Compute and print the ground state energy
  double e0 = eig0(bonds, block);
  HydraPrint(e0);
  
  return EXIT_SUCCESS;
}

```

### Documentation
The full documentation is available at [awietek.github.io/hydradoc](https://awietek.github.io/hydradoc).

### About
author:   Alexander Wietek
license:   Apache License 2.0
