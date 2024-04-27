---
title: Quick start
---

## Hello World

Let us set up our first program using the `xdiag` library. 

=== "C++"
	```c++ 
	--8<-- "examples/hello_world/main.cpp"
	```

=== "Julia"

	```julia 
	--8<-- "examples/hello_world/main.jl"
	```

The function `say_hello()` prints out a welcome message, which also contains information which exact XDiag version is used.

In Julia this is all there is to it. For the C++ code we need to create two files to compile the program. The first is the actual `C++` code. What is maybe a bit unfamiliar is the `try / catch` block. XDiag implements a traceback mechanism for runtime errors, which is activated by this idiom. While not stricly necessary here, it is a good practice to make use of this.

Now that the application program is written, we next need to set up the compilation instructions using [CMake](https://cmake.org/). To do so we create a second file called `CMakeLists.txt` in the same directory.

```cmake
--8<-- "examples/hello_world/CMakeLists.txt"
```

You should replace `"/path/to/xdiag/install"` with the appropriate directory where your XDiag library is installed after compilation. This exact `CMakeLists.txt` file can be used to compile any XDiag application.

!!! info

    For using the distributed XDiag library the last line of the above
    `CMakeLists.txt` should be changed to

    ```cmake
    target_link_libraries(main PUBLIC xdiag::xdiag_distributed)
    ```

We then compile the application code,

```bash
cmake -S . -B build
cmake --build build
```

and finally run our first `xdiag` application.

```bash
./build/main
```

## Computing the ground state energy of a spin chain 

We compute the ground state energy of the $S=1/2$ Heisenberg chain on
a periodic chain lattice in one dimension. The Hamiltonian is given by

$$ H = J\sum_{\langle i,j \rangle} \mathbf{S}_i \cdot \mathbf{S}_j$$

where $\mathbf{S}_i = (S_i^x, S_i^y, S_i^z)$ are the spin $S=1/2$ operators
and $\langle i,j \rangle$ denotes summation over nearest-meighbor sites
$i$ and $j$.

The following code, sets up the Hilbert space, defines the Hamiltonian and finally calls an iterative eigenvalue solver to compute the ground state energy.

=== "C++"
	```c++
	--8<-- "examples/spectrum/spinhalf_chain_e0/main.cpp"
	```

=== "Julia"
	```julia
	--8<-- "examples/spectrum/spinhalf_chain_e0/main.jl"
	```

