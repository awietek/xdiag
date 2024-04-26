---
title: Installation
---

## Julia Installation

Enter the package mode using `]` in the Julia REPL and add the following two packages
```julia
add https://github.com/awietek/XDiag_jll.jl.git
add https://github.com/awietek/XDiag.jl.git
```

## C++ Compilation

Using XDiag with C++ is a two-step process. First the `xdiag` library needs
to be compiled and installed. Therafter, application codes are compiled
in a second step. Here we explain how to compile the library.

### Prerequisites

* A C++ compiler that supports C++17 (`g++`, `clang`, or Intel's `icpx`)
* [git](https://git-scm.com/) version control system
* [CMake](https://cmake.org/) build system generator 
* A linear algebra backend (Blas/Lapack, IntelMKL or Accelerate on OSX)
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

### Advanced Compilation

- **Parallel compilation**
    To speed up the compilation process, the build step can be performed in parallel using the `-j` flag

    ```bash
    cmake --build build -j
    ```

- **Listing compile options**

    The available compilation options can be displayed using
    ``` bash
    cmake -L .
    ```

- **Choosing a certain compiler**

    The compiler (e.g. `icpx`) can be specified using
    ``` bash
    cmake -S . -B build -D CMAKE_CXX_COMPILER=icpx
    ```

    !!! warning 

        If the `xdiag` library is compiled with a certain compiler, it is
        advisable to also compile the application codes with the same compiler.

- **Setting the install path**

    In the installation step, the install directory can be set in the following way
    ```bash
    cmake --install build --prefix /my/install/prefix
    ```

- **Disabling HDF5/OpenMP**

    To disable support for [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
    or [OpenMP](https://www.openmp.org/) support, use
    ```bash
    cmake -S . -B build -D XDIAG_DISABLE_OPENMP=On -D XDIAG_DISABLE_HDF5=On
    ```
    
- **Building and running tests**

    To compile and run the testing programs, use
    ``` bash
    cmake -S . -B build -D BUILD_TESTING=On
    cmake --build build
    build/tests/tests
    ```

- **Building the Julia wrapper locally**

    First, get the path to the `CxxWrap` package of julia. To do so, enter the Julia REPL,
    ```bash
    julia
    ```
    and print the corresponding path using
    ```julia
    using CxxWrap
    CxxWrap.prefix_path()
    ```
    This should print the `/path/to/libcxxwrap-julia-prefix`. This is then used to configure the cmake compilation.
    ``` bash
    cmake -S . -B build -D XDIAG_JULIA_WRAPPER=On -D CMAKE_PREFIX_PATH=/path/to/libcxxwrap-julia-prefix
    cmake --build build
    cmake --install build
    ```
    The julia wrapper library can then be found in the install dir as `libxdiagjl.so`, (or the corresponding library format on non-Linux systems).
	
## Building Documentation
The source files for the documentation can be found in the directory `docs`. The documentation is built using [Material for MKDocs](https://squidfunk.github.io/mkdocs-material/). To work on it locally, it can be served using 

```bash
mkdocs serve
```

from the `xdiag` root source directory. A local build of the documentation can then be accessed in a webbrowser at the adress

```
127.0.0.1:8000
```
