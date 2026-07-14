## Installation


### Julia library

XDiag can conveniently be installed via the Julia package manager. First, a Julia interpreter needs to be opened from a command line using `julia`. Then, the "package mode" can be entered by typing `]`, and XDiag is installed using `add XDiag`. In summary:

=== "Bash"
	```bash
	$ julia
	julia> ]
	pkg> add XDiag
	```

---

### C++ library

The first step in employing the C++ version is to compile the library. The source code can be obtained from [github](https://github.com/awietek/xdiag) by cloning using [git](https://git-scm.com/).

=== "Bash"
	```bash
	cd /path/to/where/xdiag/should/be
	git clone https://github.com/awietek/xdiag.git
	```

The compilation and installation is then performed with [CMake](https://cmake.org/) 

=== "Bash"
	```bash
	cd xdiag
	cmake -S . -B build
	cmake --build build
	cmake --install build
	```

By default, the resulting library is installed at `/path/to/where/xdiag/should/be/install`. There are various options when compiling, including optimizations that can be used. For more details on the compilation process, we refer to the [compilation guide](documentation/compilation/compilation.md).
