---
title: Examples
---

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

<!-- - Time-evolution of a $t$-$J$ model -->

<!--     ??? example "source" -->
<!--         ```c++  -->
<!--         --8<-- "examples/time_evolution/tj_distributed_time_evolve/main.cpp" -->
<!--         ``` -->

<!-- ## Application CMakeLists.txt -->

<!-- - Normal xdiag library -->

<!--     ??? example "source" -->
<!--         ```cmake  -->
<!--         --8<-- "examples/hello_world/CMakeLists.txt" -->
<!--         ``` -->

<!-- - Distributed xdiag library -->

<!--     ??? example "source" -->
<!--         ```cmake  -->
<!--         --8<-- "examples/time_evolution/tj_distributed_time_evolve/CMakeLists.txt" -->
<!--         ``` -->
