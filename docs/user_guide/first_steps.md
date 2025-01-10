---
title: First steps
---

## Writing code

Let us set up our first program using the `xdiag` library. 

=== "Julia"

	```julia 
	--8<-- "examples/hello_world/main.jl"
	```
	
=== "C++"
	```c++ 
	--8<-- "examples/hello_world/main.cpp"
	```

The function `say_hello()` prints out a welcome message, which also contains information which exact XDiag version is used. In Julia this is all there is to it.

For the C++ code we need to create two files to compile the program. The first is the actual `C++` code. What is maybe a bit unfamiliar is the `try / catch` block. XDiag implements a traceback mechanism for runtime errors, which is activated by this idiom. While not stricly necessary here, it is a good practice to make use of this.

## Compilation

Now that the application program is written, we next need to set up the compilation instructions using [CMake](https://cmake.org/). To do so we create a second file called `CMakeLists.txt` in the same directory.

```cmake
--8<-- "examples/hello_world/CMakeLists.txt"
```

You should replace `"/path/to/xdiag/install"` with the appropriate directory where your XDiag library is installed after compilation. This exact `CMakeLists.txt` file can be used to compile any XDiag application.

!!! info

    For using the distributed XDiag library the last line of the above
    `CMakeLists.txt` should be changed to

    ```cmake
    target_link_libraries(main PUBLIC xdiag::xdiag_distributed)
    ```

We then compile the application code,

```bash
cmake -S . -B build
cmake --build build
```
and finally run our first XDiag application.

```bash
./build/main
```

## Hilbert spaces

We are now ready to run our first actual calculation using XDiag. Our immediate goal will be to determine the ground state energy of a $S=1/2$ Heisenberg model on a 1D chain lattice with periodic boundary conditions,

$$ H = J \sum_{\langle i, j\rangle} \mathbf{S}_i \cdot \mathbf{S}_j,$$

where $\mathbf{S}_i = (S^x_i, S^y_i, S^z_i)$ denotes the vector of spin matrices at a given site $i$. The notation $\langle i, j\rangle$ refers to summatation over neighboring sites $i$ and $j$.

The first thing to define before any computation, is the Hilbert space our model will be defined on. For the present example, we use the Hilbert space class [Spinhalf](../documentation/blocks/spinhalf.md). We consider a chain lattice with $N=8$ sites and create a $S=1/2$ Hilbert space:

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:first_steps_1"
	```

We would like to know which spin configurations, the Hilbert space is made up of. To do so, we can iterate over the Hilbert space and print out the spin configurations. 


=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:first_steps_2"
	```

This produces an output similar to the following:

```bash
↓↓↓↓↓↓↓↓
↓↓↓↓↓↓↓↑
↓↓↓↓↓↓↑↓
↓↓↓↓↓↓↑↑
↓↓↓↓↓↑↓↓
↓↓↓↓↓↑↓↑
↓↓↓↓↓↑↑↓
...
```

Here we already see several things at work. XDiag features a convenient way to write logs, with the [Log](../documentation/utilities/logging.md) class. The first argument to `Log()` is a format string. In C++ we use the [fmt](https://fmt.dev/) library, to be able to write structured format and format our output. The second argument turns our `spins` into a string. `spins` is of type [ProductState](../documentation/states/product_state.md), whose configuration on each site can be individually addressed. 

Further, we notice that all $2^N$ spin configurations are included in this Hilbert space. However, the Heisenberg model conserves the total $S^z = \sum_i S^z_i$, and thus we could limit ourselves to a block of the Hilbert space, which only contains configurations of a certain magnetization: 
	
=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:first_steps_3"
	```
	
This produces an output similar to:
```bash
↓↓↓↓↑↑↑↑
↓↓↓↑↓↑↑↑
↓↓↓↑↑↓↑↑
↓↓↓↑↑↑↓↑
↓↓↓↑↑↑↑↓
↓↓↑↓↓↑↑↑
↓↓↑↓↑↓↑↑
...
```
	
We see that the block now only contains configurations, where the total number of spins pointing up is 4. We can now quickly check the dimension of the Hilbert spaces, and confirm that the dimension of the block is reduced from $2^8=256$ to $\begin{pmatrix} 8 \\ 4 \end{pmatrix} = 70$,

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:first_steps_4"
	```

which should print:
```bash
hspace.size():
256
block.size():
70
```

Here, we introduced another functionality of XDiag, the [XDIAG_SHOW](../documentation/utilities/xdiag_show.md) macro which can be used for quick debug printing of XDiag objects.
