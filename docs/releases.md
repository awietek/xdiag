---
title: Releases
---

## v0.3.1

**Feb. 11, 2025**

New API and operator logic

* Updated to a new streamlined API
* Implemented operator logic including symmetries
* Completed wrapper for Julia
* w = apply(ops, v) works now

---

## v0.2.3

**Sep. 9, 2024**

Introduced 1-indexing everywhere in Julia version

* only changes to XDiag.jl, C++ untouched
* XDiag_jll.jl remains at v0.2.2

---

## v0.2.2

**Aug. 27, 2024**

Lanczos routines and multicolumn States

* wrapped eigs_lanczos, eigvals_lanczos
* implemented apply for States with multiple columns
* changed wrapping of std::vectors of Op and Permutation

---

## v0.2.1

**Aug. 16, 2024**

Small patch release providing small utility functions

* wrapped say_hello, print_version, and set_verbosity
* resorted to compiling wrapper with conventional OpenMP on aarch64 apple
* Fixed faulty behaviour of OpenMP on aarch64 apple

---

## v0.2.0

**Aug. 15, 2024**

Basic functionality for three Hilbert space types, Spinhalf, tJ, and Electron, has been implemented. Features are:

* Algebra with and without permutation symmetries
* Parallelization with OpenMP and MPI
* CMake has been properly set up
* Iterative algorithms present, Lanczos, Arnoldi, time evolution
* A minimal Julia wrapper has been written
* The Julia wrapper compiles on several target 64bit systems using BinaryBuilder
