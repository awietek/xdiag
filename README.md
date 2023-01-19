![license](https://img.shields.io/github/license/awietek/hydra)
![cpp](https://img.shields.io/badge/C++-17-blue.svg)
[![docs](https://img.shields.io/badge/Documentation-here-red.svg)](https://awietek.github.io/hydradoc)
[![Linux CI](https://github.com/awietek/hydra/actions/workflows/linux.yml/badge.svg?style=for-the-badge)](https://github.com/awietek/hydra/actions/workflows/linux.yml)
[![Mac OSX CI](https://github.com/awietek/hydra/actions/workflows/macosx.yml/badge.svg?style=for-the-badge)](https://github.com/awietek/hydra/actions/workflows/macosx.yml)
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

### Example Code:
```cpp
// Hubbard model on a N=8 linear chain with U/t = 5
int n_sites = 8;
int nup = 4;
int ndn = 4;

// Create the model, its space group and irreps
auto [bondlist, couplings] = get_linear_chain(n_sites, 1.0, 5.0);
auto [space_group, irreps] = get_cyclic_group_irreps(n_sites);

// Compute the Hamiltonian matrix
auto block = Electron(n_sites, nup, ndn, space_group, irrep);
auto H = MatrixReal(bondlist, couplings, block, block);

// Compute its eigenvalues and print them
auto evals = lila::EigenvaluesSym(H);
LilaPrint(evals);
```

### Documentation
The full documentation is available at [awietek.github.io/hydradoc](https://awietek.github.io/hydradoc).

### About
author:   Alexander Wietek
license:   Apache License 2.0
