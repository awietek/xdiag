---
title: Installation
---

## Julia installation

Enter the package mode using `]` in the Julia REPL and type:
```julia
add XDiag
```

That's it!

## C++ compilation

Using XDiag with C++ is a two-step process. First the `xdiag` library needs
to be compiled and installed. Therafter, application codes are compiled
in a second step. The library can be compiled in two different versions:

- **Normal** library: features parallelization using OpenMP only
- **Distributed** library: features distributed-memory parallelization using MPI. 

### Prerequisites

* A C++ compiler that supports C++17 (`g++`, `clang`, or Intel's `icpx`)
* [git](https://git-scm.com/) version control system
* [CMake](https://cmake.org/) build system generator 
* A linear algebra backend (BLAS/LAPACK, Intel MKL or Accelerate on OSX)
* **optional** [HDF5](https://www.hdfgroup.org/solutions/hdf5/), [OpenMP](https://www.openmp.org/)
* **optional** [MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface) for the distributed library


### Basic Compilation

- **Download** the source code using [git](https://git-scm.com/)
  ```bash
  cd /path/to/where/xdiag/should/be
  git clone https://github.com/awietek/xdiag.git
  ```

- **Compile the default library**
  ``` bash
  cd xdiag
  cmake -S . -B build
  cmake --build build
  cmake --install build
  ```
  By default, the library is now installed in the subdirectory `install`.

- **Compile the distributed library**

    To use the distributed computing features of `xdiag`, the distributed
    library has to be built which requires [MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface).
    ``` bash
    cd xdiag
    cmake -S . -B build -D XDIAG_DISTRIBUTED=On
    cmake --build build
    cmake --install build
    ```

    !!! info

        It might be necessary to explicitly define MPI compiler, e.g. `mpicxx` like this
        ```bash
        cmake -S . -B build -D XDIAG_DISTRIBUTED=On -D CMAKE_CXX_COMPILER=mpicxx
        ```
		
		
The compilation process can be modified and also allows for further optimizations. We collect several common scenarios in the [Advanced Compilation](../documentation/compilation/advanced_compilation.md) guide.
