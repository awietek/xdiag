---
title: Examples
---

## Basic examples

- Hello World

    ??? example "source"
        ```c++ 
        --8<-- "examples/hello_world/main.cpp"
        ```
	
- Ground state energy of a spin chain

    ??? example "source"
        ```c++ 
        --8<-- "examples/spectrum/spinhalf_chain_e0/main.cpp"
        ```

## Distributed examples


- Time-evolution of a $t$-$J$ model

    ??? example "source"
        ```c++ 
        --8<-- "examples/time_evolution/tj_distributed_time_evolve/main.cpp"
        ```

## Application CMakeLists.txt

- Normal hydra library

    ??? example "source"
        ```cmake 
        --8<-- "examples/hello_world/CMakeLists.txt"
        ```

- Distributed hydra library

    ??? example "source"
        ```cmake 
        --8<-- "examples/time_evolution/tj_distributed_time_evolve/CMakeLists.txt"
        ```