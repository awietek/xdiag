### Writing application code

To employ the XDiag library, an application code is written. The simplest application using XDiag is a `hello_world` program.

=== "C++"
	```c++ 
	--8<-- "examples/hello_world/main.cpp"
	```
=== "Julia"
	```julia 
	--8<-- "examples/hello_world/main.jl"
	```

The function `say_hello()` prints out a welcome message, which also contains information on which exact XDiag version is used. We would like to emphasize the `try / catch` block in the C++ version. XDiag implements a traceback mechanism for runtime errors, which is activated by the `error_trace` function. While not strictly necessary here, it is a good practice to employ this.

---

### Application compilation

In C++, now that the application program is written, we next need to set up the compilation instructions using [CMake](https://cmake.org/). To do so we create a second file called `CMakeLists.txt` in the same directory.

=== "CMake"
	```cmake
	--8<-- "examples/hello_world/CMakeLists.txt"
	```

You should replace `/path/to/where/xdiag/should/be/install` with the appropriate directory where your XDiag library is installed after compilation. We then compile the application code,

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
