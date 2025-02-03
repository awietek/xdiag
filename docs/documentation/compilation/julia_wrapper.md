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


In order to test the new shared library `xdiagjl` together with the `XDiag.jl` library, we need to override the artifact associated with the `XDiag_jll.jl` package. For this, first we need to find out which artifact is associated with `XDiag_jll.jl`. For this enter julia and type the following commands:

``` bash
julia
using XDiag_jll
XDiag_jll.find_artifact_dir()
```

This gives the directory in which the artifact is defined. Now we have to add a line in the `Overrides.toml` file, typically located at `.julia/artifacts/Overrides.toml`.

Here, we then add a line like this:

```toml
55ec928f6054024a4e9bf02e74e4da8b69175655 = "/path/to/xdiag/install"
```

The hash is to be replaced by the directory of the `XDiag_jll.jl` artifact.
