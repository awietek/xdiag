---
title: Releases
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
