---
title: Installation
---

## Julia
### Installing latest release
```julia
add YxDiag
```

### Compiling the julia shared library
```julia
add YxDiag
```

## C++

### Prerequisites
* A C++ compiler that supports C++17, e.g. GNU's `g++`, `clang`, or Intel's `icpx`
* [CMake](https://cmake.org/) build system generator 
* A linear algebra backend, e.g. the Netlib Blas/Lapack, IntelMKL or Accelerate on OSX.
* [git](https://git-scm.com/) version control system


### Downloading the source code

**Step 1:**
Choose an installation directory

```bash
cd /path/to/hydra
```

**Step 2:**
Clone the hydra library from github

```bash
git clone https://github.com/awietek/hydra.git
cd hydra
```

### Compiling the default library

``` bash
cd /path/to/hydra/dir
cmake -S . -B build
cmake --build build
cmake --install
```

### Compiling the distributed library

``` bash
cd /path/to/hydra/dir
cmake -S . -B build -D BUILD_DISTRIBUTED=On
cmake --build build
cmake --install
```