### Distributed computing

The standard XDiag library features shared memory parallelization using [OpenMP](https://www.openmp.org/) for both the C++ and Julia version. However, we also provide distributed memory parallelization in a separate C++-only library that needs to be compiled independently. To do so, we have to use the flag `XDIAG_DISTRIBUTED` when setting up the compilation using [CMake](https://cmake.org/).

```bash
cmake -S . -B build -D XDIAG_DISTRIBUTED=On
cmake --build build
cmake --install build
```

This will create a distinct `xdiag_distributed` library, which is different from the standard `xdiag` library. To link an application code to the distributed library, we can use the following `CMakeLists.txt` file.

```cmake
cmake_minimum_required(VERSION 3.19)
project(hello_world)
find_package(xdiag_distributed REQUIRED HINTS "/path/to/where/xdiag/should/be/install")
add_executable(main main.cpp)
target_link_libraries(main PRIVATE xdiag::xdiag_distributed)  
```

Notice that this only differs from the original `CMakeLists.txt` file shown in the section [Application compilation](#application-compilation) in two places, where we specify that the `xdiag_distributed` library instead of the standard `xdiag` library is used. The distributed memory library is built on top of the message passing interface (MPI). Every MPI application needs to initialize and finalize the MPI environment explicitly. Hence, a typical main routine for using the distributed XDiag library with MPI should look similar to the following listing.

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:usage_guide_dist1"
	```

The functionality described in the previous sections is also available to some degree for the distributed library. Importantly, to use the distributed capabilities we have to change the types of blocks from the standard [Spinhalf](documentation/blocks/spinhalf.md), [tJ](documentation/blocks/tJ.md), and [Electron](documentation/blocks/electron.md) blocks to the [SpinhalfDistributed](documentation/blocks/spinhalf_distributed.md), [tJDistributed](documentation/blocks/tJ_distributed.md), and [ElectronDistributed](documentation/blocks/electron_distributed.md) blocks. An example to compute the ground state of a Heisenberg chain using the distributed capabilities is given by,

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:usage_guide_dist2"
	```

At present, one important difference between the distributed and standard library is that the distributed blocks do not yet support symmetrized blocks. This feature will be added in future versions of XDiag.
