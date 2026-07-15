---
title: Installation
---

# Installation

XDiag comes in two flavors: a [Julia](#julia-library) package for interactive,
high-level usage and a [C++](#c-library) library for maximal performance and
access to the distributed-memory features. Both share the same API, so the
code examples throughout this user guide are given for either language.

## Julia library

XDiag can conveniently be installed via the Julia package manager. First, a Julia interpreter needs to be opened from a command line using `julia`. Then, the "package mode" can be entered by typing `]`, and XDiag is installed using `add XDiag`. In summary:

=== "Julia"
	```julia
	julia> ]
	pkg> add XDiag
	```

That's it! No further compilation step is required.

## C++ library

The first step in employing the C++ version is to compile the library. The source code can be obtained from [GitHub](https://github.com/awietek/xdiag) by cloning it using [git](https://git-scm.com/).

=== "Bash"
	```bash
	cd /path/to/where/xdiag/should/be
	git clone https://github.com/awietek/xdiag.git
	```

The compilation and installation is then performed with [CMake](https://cmake.org/).

=== "Bash"
	```bash
	cd xdiag
	cmake -S . -B build
	cmake --build build
	cmake --install build
	```

By default, the resulting library is installed in the `install` subdirectory of the source tree. There are various options when compiling, including performance optimizations, HDF5/OpenMP support, and building the distributed library with MPI.

!!! info

	For the full list of prerequisites, compile-time options, the distributed
	build, and how to compile your application code against XDiag, see the
	[Installation](../installation.md) page and the detailed
	[compilation guide](../documentation/compilation/compilation.md).
