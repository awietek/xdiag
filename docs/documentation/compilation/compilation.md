---
title: Compilation
---

## Library Compilation

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

### Download

The source code can be downloaded from Github using [git](https://git-scm.com/) 
```bash
cd /path/to/where/xdiag/should/be
git clone https://github.com/awietek/xdiag.git
```

### Default library
``` bash
cd xdiag
cmake -S . -B build
cmake --build build
cmake --install build
```
By default, the library is now installed in the subdirectory `install` of your XDiag path `/path/to/where/xdiag/should/be`. If you would like to install it to another directory, you could set up the CMake compilation using,
  
```bash
cmake -S . -B build -D CMAKE_INSTALL_PREFIX=/path/to/where/you/want/xdiag/installed
cmake --build build
cmake --install build
```

### Distributed library

To use the distributed computing features of `xdiag`, the distributed library has to be built which requires [MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface).
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
		
## Application Compilation

Once an application program is written, we next need to set up the compilation instructions using [CMake](https://cmake.org/). To do so we a second file called `CMakeLists.txt` in the application directory.

```cmake
--8<-- "examples/cmake/CMakeLists.txt"
```

You should replace `"/path/to/xdiag/install"` with the appropriate directory where your XDiag library is installed after compilation. This exact `CMakeLists.txt` file can be used with little alteration to compile any XDiag application.

!!! info

    For using the distributed XDiag library the above `CMakeLists.txt` should be changed to.

    ```cmake
    cmake_minimum_required(VERSION 3.19)

	project(xdiag_application)

	find_package(xdiag_distributed REQUIRED HINTS "/path/to/xdiag/install")
	add_executable(main main.cpp)
	target_link_libraries(main PRIVATE xdiag::xdiag_distributed)
    ```
	
	Notice that we only need to replace `xdiag` by `xdiag_distributed` in two places.

We then compile the application code,

```bash
cmake -S . -B build
cmake --build build
```

and finally run our XDiag application.

```bash
./build/main
```
		
## Advanced Compilation

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

- **Shared / Static libraries**

    XDiag can be built either as a static or shared library. By default, a static library is
    built. To build a shared library, use the option **XDIAG_SHARED_LIBS**, e.g.
    ``` bash
    cmake -S . -B build -D XDIAG_SHARED_LIBS=On
    ```
    Typically, shared libraries can reduce the size of the executable and are often prefered.
    However, when building an application and linking to a static library, the resulting
    binary can be used standalone, without the shared library needing to be loaded. This means
    once an application code is compiled, it will run the same way indefinitely, even if the
    XDiag library changes. Thus, it is mostly more convenient to link against the static library.

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

- **Building examples**

    To compile and run the example programs, use
    ``` bash
    cmake -S . -B build -D BUILD_EXAMPLES=On
    cmake --build build
    build/examples/usage_examples
    ```


## Optimization


- **Native optimization**

	Adds the flag `-march=native` to perform optimizations for the native architecture. This can have tremendous performance impact, especially on t-J models, since then the BMI2 instructions are activated whenever they are available (e.g. newer Intel and AMD processors)
    ``` bash
    cmake -S . -B build -D XDIAG_OPTIMIZE_FOR_NATIVE=On
    ```
