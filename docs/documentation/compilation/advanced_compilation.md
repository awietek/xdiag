---
title: Advanced Compilation
---


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
