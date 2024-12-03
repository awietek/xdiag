---
title: Building the Julia wrapper
---

In order to develop and extend the julia wrapper, one should work locally and build a local version of the `xdiag` Julia binaries. First, get the path to the `CxxWrap` package of julia. To do so, enter the Julia REPL,
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
