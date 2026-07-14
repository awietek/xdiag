<style>
  .md-typeset h1,
  .md-content__button {
    display: none;
  }
</style>

<img src="img/logo_cropped.png" alt="XDiag logo" width="500"/>

**Exact Diagonalization of quantum many-body systems — in Julia and C++.**

[Quick Start](quick_start.md){ .md-button .md-button--primary }
[User Guide](user_guide/index.md){ .md-button .md-button--secondary }
[Publications using XDiag](publications.md){ .md-button .md-button--secondary }

| Languages                                                                                                                                                                                                      | Code status                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |               Version |
|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------:|
| ![C++ Badge](https://img.shields.io/badge/C%2B%2B-00599C?logo=cplusplus&logoColor=fff&style=for-the-badge) ![Julia](https://img.shields.io/badge/-Julia-9558B2?style=for-the-badge&logo=julia&logoColor=white) | [![SciPost](https://img.shields.io/badge/Publication-SciPostPhysCodeb.70-yellow)](https://scipost.org/10.21468/SciPostPhysCodeb.70)[![Linux CI](https://github.com/awietek/xdiag/actions/workflows/linux.yml/badge.svg?style=for-the-badge)](https://github.com/awietek/xdiag/actions/workflows/linux.yml) [![Mac OSX CI](https://github.com/awietek/xdiag/actions/workflows/osx.yml/badge.svg?style=for-the-badge)](https://github.com/awietek/xdiag/actions/workflows/osx.yml)[![Intel MPI CI](https://github.com/awietek/xdiag/actions/workflows/intelmpi.yml/badge.svg?style=for-the-badge)](https://github.com/awietek/xdiag/actions/workflows/intelmpi.yml)[![Examples CI](https://github.com/awietek/xdiag/actions/workflows/examples.yml/badge.svg?style=for-the-badge)](https://github.com/awietek/xdiag/actions/workflows/examples.yml)[![Julia CI](https://github.com/awietek/XDiag.jl/actions/workflows/CI.yml/badge.svg?style=for-the-badge)](https://github.com/awietek/XDiag.jl/actions/workflows/CI.yml)  [![codecov](https://codecov.io/gh/awietek/xdiag/graph/badge.svg)](https://codecov.io/gh/awietek/xdiag) ![license](https://img.shields.io/badge/license-Apache%202.0-blue) | [v0.5.0](releases.md) |

## Overview

XDiag is a library for performing **exact diagonalizations** of quantum many-body systems. It is designed to be both easy to use and highly performant, and consists of two packages that share a common API:

* the core C++ library [xdiag](https://github.com/awietek/xdiag), and
* the convenient Julia wrapper [XDiag.jl](https://github.com/awietek/XDiag.jl).

Key features:

* **Many Hilbert space types** — spins, electrons, $t$-$J$ models, spinless fermions, and bosons.
* **A full operator algebra** — assemble arbitrary operators and multiply or hermitian-conjugate them freely.
* **Iterative algorithms** — Lanczos and LOBPCG diagonalization, and real- or imaginary-time evolution.
* **Symmetries** — exploit space-group symmetries through symmetry-adapted bases.
* **Parallelization** — shared-memory with OpenMP and distributed-memory with MPI.
* **Optimized combinatorics** for fast navigation of large Hilbert spaces.

## News
* **Jul. 30, 2026** — New **major** release [v0.5.0](releases.md): a generic operator algebra, spinless-fermion and boson blocks, an LOBPCG eigensolver, and improved memory efficiency.
* **Dec. 4, 2025** — New release [v0.4.1](releases.md) fixes a bug reported as issue #99.
* **Nov. 4, 2025** — New release [v0.4.0](releases.md) introduces sparse matrix capabilities.
* **Jun. 17, 2025** — New release [v0.3.3](releases.md) features additional parallelization in Julia, compatibility with Julia 1.12, and several practical enhancements.

## Gallery

A few results produced with XDiag:

<div class="grid cards" markdown>
- ![Triangular-lattice moiré model of WSe₂](img/triangularwse2.png){ align=left }
- ![Energy level statistics of a spin-1/2 chain](img/spinhalf_chain_level_statistics.png){ align=left }
- ![Green's function of the Hubbard model](img/hubbard_greens_f.png){ align=left }
- ![Excitation spectra of the J₁–J₂ model](img/j1j2_spectra.png){ align=left }
</div>
