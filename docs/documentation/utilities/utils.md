---
title: Miscellaneous utility functions
---

## set_verbosity

Set how much logging is generated my XDiag to monitor the progress and behaviour of the code. There are three verbosity levels that can be set:

* 0: no output at all
* 1: some output 
* 2: detailed output 

This can be useful, e.g. to monitor the progress of an iterative algorithm

=== "Julia"
	```julia
	set_verbosity(level::Int64)
	```

=== "C++"
	```c++
    void set_verbosity(int64_t level);
	```


## say_hello

Prints a nice welcome message containing the version number and git commit used.

=== "Julia"
	```julia
	say_hello()
	```

=== "C++"
	```c++
    void say_hello()
	```

## print_version

If `say_hello` is too much flower power for you, one can also just have a boring print-out of the version number using this function. 

=== "Julia"
	```julia
	print_version()
	```

=== "C++"
	```c++
    void print_version()
	```
