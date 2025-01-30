---
title: Compilation
---

## Basic Compilation

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

## Optimization

**missing documentation**
