---
title: Writing code
---

# Writing code

Once XDiag is [installed](01-installation.md), we are ready to write our first
application. The simplest possible XDiag program is a `hello_world` program,
which prints out a welcome message together with information on the exact XDiag
version being used. Below we describe how to set up such a program, both as a
local Julia project and as a compiled C++ application.

## Setting up a Julia project

It is good practice to give every XDiag computation its own Julia *project*
(also called an *environment*). A project keeps track of the exact package
versions used, which makes results reproducible and avoids interference between
different computations. A project is defined by a `Project.toml` and a
`Manifest.toml` file, which Julia creates and manages automatically.

We first create a directory for our project and enter it,

=== "Bash"
	```bash
	mkdir my_xdiag_project
	cd my_xdiag_project
	```

We then start a Julia interpreter with this directory activated as the current
project using the `--project` flag,

=== "Bash"
	```bash
	julia --project=.
	```

Now we enter the package mode by typing `]` and add XDiag to the project,

=== "Julia"
	```julia
	pkg> add XDiag
	```

This downloads XDiag and records it in the `Project.toml` and `Manifest.toml`
files inside our directory. We only need to do this once per project. We can now
leave the package mode by pressing backspace and write our first line of XDiag
code directly in the REPL,

=== "Julia"
	```julia
	--8<-- "examples/hello_world/main.jl"
	```

More commonly, we write our code to a file, say `main.jl`, and execute it from
the command line with the project activated,

=== "Bash"
	```bash
	julia --project=. main.jl
	```

The function `say_hello()` prints out a welcome message, which also contains
information on which exact XDiag version is used.

---

## Writing a C++ application

To employ the C++ library, an application code is written in a `main.cpp` file.
The `hello_world` program looks as follows.

=== "C++"
	```c++
	--8<-- "examples/hello_world/main.cpp"
	```

As in the Julia version, the function `say_hello()` prints out a welcome message
containing the XDiag version. We would like to emphasize the `try / catch` block
in the C++ version. XDiag implements a traceback mechanism for runtime errors,
which is activated by the `error_trace` function. While not strictly necessary
here, it is good practice to employ this in every XDiag application.

### Application compilation

Now that the application program is written, we next need to set up the
compilation instructions using [CMake](https://cmake.org/). To do so we create a
second file called `CMakeLists.txt` in the same directory.

=== "CMake"
	```cmake
	--8<-- "examples/hello_world/CMakeLists.txt"
	```

You should replace `/path/to/where/xdiag/should/be/install` with the appropriate
directory where your XDiag library is installed after compilation. We then
compile the application code,

=== "Bash"
	```bash
	cmake -S . -B build
	cmake --build build
	```

and finally run our first XDiag application.

=== "Bash"
	```bash
	./build/main
	```

!!! info

	For more details on compiling application code, including linking against the
	distributed library, see the [compilation guide](../documentation/compilation/compilation.md#application-compilation).
