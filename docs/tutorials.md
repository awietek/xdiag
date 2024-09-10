---
title: Tutorials
---

## Introduction to Exact Diagonalization using XDiag

Supporting material for lecture held at [Quant24](https://www.pks.mpg.de/de/quant24) master's school at MPI PKS. Consists of a **Jupyter notebook** and a sample lattice file describing the $N=12$ site triangular lattice Heisenberg model:

- [ed_basic_tutorial.ipynb](examples/ed_basic_tutorial.ipynb){:download="ed_basic_tutorial.ipynb"}
- [triangular.12.J1J2.toml](examples/triangular.12.J1J2.toml){:download="triangular.12.J1J2.toml"}

This notebook uses the Julia verision of XDiag and covers the basic functionality:

- How to define a Hilbert space
- How to define an operator
- How to perform a full diagonalization
- How to use $S^z$ conservation
- How to use translational symmetry
- How to use iterative algorithms for sparse diagonalization
- How to compute ground state observables

## Basic examples

<div class="grid cards" markdown>

-   :material-file-document:{ .lg .middle } __Hello World__

    ---

    Prints out a greeting containing information on the version of the code.

    [source](examples/hello_world.md) :simple-cplusplus: :simple-julia:

-   :material-file-document:{ .lg .middle } __Groundstate energy__

    ---

    Computes the ground state energy of a simple Heisenberg spin $S=1/2$ chain

    [source](examples/spinhalf_chain_e0.md) :simple-cplusplus: :simple-julia:

</div>

## Distributed examples

<div class="grid cards" markdown>

-   :material-file-document:{ .lg .middle } __$t$-$J$ time evolution__

    ---

    Computes the time evolution of a state in the $t$-$J$ model with distributed parallelization

    [source](examples/tj_distributed_time_evolve.md) :simple-cplusplus: 

</div>

## CMakeLists.txt for applications

<div class="grid cards" markdown>

-   :material-file-document:{ .lg .middle } __Normal XDiag__

    ---

    Template `CMakeLists.txt` which can be used to compile applications with the **normal** XDiag library.

    [source](examples/cmake_normal.md) :simple-cmake: 

-   :material-file-document:{ .lg .middle } __Distributed XDiag__

    ---

    Template `CMakeLists.txt` which can be used to compile applications with the **distributed** XDiag library.
	
    [source](examples/cmake_distributed.md) :simple-cmake:

</div>
