---
title: User guide
---

**A step-by-step guide to using XDiag**


## Installation


### Julia

Enter the package mode using `]` in the Julia REPL and type:
```julia
add XDiag
```

That's it!

---

### C++

For using the C++ version, we need to compile the library first. To do so, we first download the code using [git](https://git-scm.com/), 

```bash
cd /path/to/where/xdiag/should/be
git clone https://github.com/awietek/xdiag.git
```

and then compile the library using [CMake](https://cmake.org/),

``` bash
cd xdiag
cmake -S . -B build
cmake --build build
cmake --install build
```

The resulting library is now installed at `/path/to/where/xdiag/should/be/install`. There are various options when compiling, including optimizations which can be used. For more details on the compilation process, we refer to the [Compilation](documentation/compilation/compilation.md) guide.

---

## First steps

### Writing code

Let us set up our first program using the `xdiag` library. 
	
=== "C++"
	```c++ 
	--8<-- "examples/hello_world/main.cpp"
	```
=== "Julia"

	```julia 
	--8<-- "examples/hello_world/main.jl"
	```
	
The function `say_hello()` prints out a welcome message, which also contains information which exact XDiag version is used. In Julia this is all there is to it.

What is maybe a bit unfamiliar is the `try / catch` block in C++. XDiag implements a traceback mechanism for runtime errors, which is activated by the `error_trace` function. While not stricly necessary here, it is a good practice to make use of this.

---

### Compilation

In C++, now that the application program is written, we next need to set up the compilation instructions using [CMake](https://cmake.org/). To do so we create a second file called `CMakeLists.txt` in the same directory.

```cmake
--8<-- "examples/hello_world/CMakeLists.txt"
```

You should replace `/path/to/where/xdiag/should/be/install` with the appropriate directory where your XDiag library is installed after compilation. We then compile the application code,

```bash
cmake -S . -B build
cmake --build build
```
and finally run our first XDiag application.

```bash
./build/main
```

---

### Hilbert spaces

We are now ready to run our first actual calculation using XDiag. Our immediate goal will be to determine the ground state energy of a $S=1/2$ Heisenberg model on a 1D chain lattice with periodic boundary conditions,

$$ H = J \sum_{\langle i, j\rangle} \mathbf{S}_i \cdot \mathbf{S}_j,$$

where $\mathbf{S}_i = (S^x_i, S^y_i, S^z_i)$ denotes the vector of spin matrices at a given site $i$. The notation $\langle i, j\rangle$ refers to summatation over neighboring sites $i$ and $j$.

The first thing to define before any computation, is the Hilbert space our model will be defined on. For the present example, we use the Hilbert space class [Spinhalf](documentation/blocks/spinhalf.md). Further possible Hilbert spaces include [tJ](documentation/blocks/tJ.md) and [Electron](documentation/blocks/electron.md), see [Blocks](documentation/index.md#blocks). We consider a chain lattice with $N=8$ sites and create a $S=1/2$ Hilbert space:

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:first_steps_1"
	```
=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:first_steps_1"
	```
	
We would like to know which spin configurations, the Hilbert space is made up of. To do so, we can iterate over the Hilbert space and print out the spin configurations. 


=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:first_steps_2"
	```
=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:first_steps_2"
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

Here we already see several things at work. XDiag features a convenient way to write logs in C++, with the [Log](documentation/utilities/logging.md) class. The first argument to `Log()` is a format string. In C++ we use the [fmt](https://fmt.dev/) library, to be able to write structured format and format our output. The second argument turns our `spins` into a string. `spins` is of type [ProductState](documentation/states/product_state.md), whose configuration on each site can be individually addressed. 

Further, we notice that all $2^N$ spin configurations are included in this Hilbert space. However, the Heisenberg model conserves the total $S^z = \sum_i S^z_i$, and thus we could limit ourselves to a block of the Hilbert space, which only contains configurations of a certain magnetization: 
	
	
=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:first_steps_3"
	```
=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:first_steps_3"
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
=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:first_steps_4"
	```
	
which should print:
```bash
size(hspace):
256
size(block):
70
```

Here, we introduced another functionality of XDiag in C++, the [XDIAG_SHOW](documentation/utilities/xdiag_show.md) macro which can be used for quick debug printing of XDiag objects.

---

### Operators

Next, we define our Hamiltonian. We do so by using an [OpSum](documentation/operators/opsum.md) object. These encode sums of operators of the form 

$$ \mathcal{O} = \sum_i c_i \mathcal{O}_i. $$

Here, $\mathcal{O}_i$ denote single operators described as [Op](documentation/operators/op.md) objects and $c_i$ are the coupling constants. These can either be a real or complex number, or a string which later needs to be replaced. THe Hamiltonian of the spin $S=1/2$ Heisenberg chain is created like this:

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:first_steps_5"
	```
	
=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:first_steps_5"
	```

We first create an empty [OpSum](documentation/operators/opsum.md) and then add additional terms to it. The first part of the product denotes the coupling constant, here given as a string. Alternatively, one could have directly used real / complex numbers here. The second part of the product is a single [Op](documentation/operators/op.md) object. It is created with two inputs:

1. The type, here `SdotS` denoting an operator of the form $\mathbf{S}_i\cdot\mathbf{S}_j$. XDiag features a wide variety of operator types, see [Operator types](documentation/operators/operator_types.md).
2. An array defining which site the operator lives on. Notice, that in julia we start counting the sites from 1, while in C++ we start counting the sites from 0.

---

### Ground states

Now our Hamiltonian is defined, we can directly compute the ground state and its energy using the function [eig0](documentation/algorithms/eig0.md):

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:first_steps_6"
	```
=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:first_steps_6"
	```
	
Here, `e0` is a double precision real number and `psi0` is a [State](documentation/states/state.md) object. The [eig0](documentation/algorithms/eig0.md) uses an iterative algorithms to compute the ground state by applying the Hamiltonian in a matrix-free manner, and is thus very useful for large system sizes. It is based on a Lanczos algorithm, which can also be called directly using [eigs_lanczos](documentation/algorithms/eigs_lanczos.md) which offers more control on the behavior of the algorithm.

Alternatively, we can also perform a full diagonalization by computing the full Hamiltonian matrix using the [matrix](documentation/algebra/matrix.md) function.

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:first_steps_7"
	```
=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:first_steps_7"
	```

Notice, that we are using the [Armadillo](https://arma.sourceforge.net) library with the `arma` namespace in C++. This library serves as the linear algebra backend of XDiag and can be used to perform further calculations. In Julia, the `eigen` and `Symmetric` functions are part of the `LinearAlgebra` standard library.

---

### Measurements

Finally, we might be interested in measuring some observables of the ground state. We compute the static spin correlation $S^z_0 S^z_j$:

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:first_steps_8"
	```
=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:first_steps_8"
	```
	
Here, we are using the [inner](documentation/algebra/algebra.md#inner) function to compute an expectation value of the form $\langle \psi_0 |\mathcal{O}| \psi_0 \rangle$ with $\mathcal{O} = S^z_0 S^z_j$.

---

## Input / Output

Julia features a variety of packages facilitating input and output of data. For C++, XDiag provides convenient functionality for [TOML](https://toml.io/en/) and [HDF5](https://www.hdfgroup.org/solutions/hdf5/) files. 


### Reading from TOML

While defining a Hamiltonian or other operators in code as above can be done, it is often preferable to define operators in a file and read it in. In XDiag, we use the TOML language to define input. The Hamiltonian of the $N=8$ site Heisenberg chain we created above can be written in a TOML file as a list like this:

```toml
--8<-- "examples/user_guide/spinhalf_chain.toml:ops"
```

The first entry in every list element is the coupling constant `J`, the second enty is the type `SdotS`, and the following two entries are the sites of the operator. To read in such a Hamiltonian from a toml file we can use the [FileToml](documentation/io/file_toml.md) object together with the [read_opsum](documentation/io/read_opsum.md) Function,

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:io_1"
	```
=== "Julia"
	```julia 
	fl = FileToml("spinhalf_chain.toml")
	ops = read_opsum(fl, "Interactions")
	```
	
XDiag also features the functions [read_permutation_group](documentation/io/read_permutation_group.md) and [read_representation](documentation/io/read_representation.md) to conventiently read [PermutationGroup](documentation/symmetries/permutation_group.md) and [Representation](documentation/symmetries/representation.md) objects used to describe symmetries from file, see further below.
	
---

### Writing results to hdf5

After a simulation, we want to store results to file. A standard scientific data format is the [hdf5](https://www.hdfgroup.org/solutions/hdf5/) format. Julia supports input and output to hdf5 with the [HDF5.jl](https://juliaio.github.io/HDF5.jl/stable/) package. 

For C++ we provide a convenient way of writing results to hdf5 files. In general all numerical data, including scalar real/complex numbers as well as armadillo vectors and matrices can be easily written to file using the [FileH5](documentation/io/file_h5.md) object.

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:io_2"
	```

The second argument `"w!"` specifies the access mode to the file. `"w!"` is the forced write mode, where preexisting files are overwritten. 

---

## Symmetries

### Permutations

Key functionality of XDiag comprises symmetry-adapted calculations, in particular with permutation symmtries which can be translation symmetries on a lattice. The Heisenberg chain with periodic boundary conditions is invariant under the translation operator

$$ T: i \rightarrow i+1 $$

To represent such symmetries, we use the [Permutation](documentation/symmetries/permutation.md) object. For example, on the $8$-site chain lattice, the symmetry which maps the sites, 

$$ \{0,1,2,3,4,5,6,7\} \rightarrow \{1,2,3,4,5,6,7,0\} $$

can be represented via: 

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:symmetries_1"
	```
=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:symmetries_1"
	```

Notice, that also here we start counting from 1 in Julia, and from 0 in C++. [Permutation](documentation/symmetries/permutation.md) objects can be multiplied, inverted, and raised to a power. We can use them to define a [PermutationGroup](documentation/symmetries/permutation_group.md) object in the following way:

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:symmetries_2"
	```
=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:symmetries_2"
	```

### Representations

One-dimensional representations of a [PermutationGroup](documentation/symmetries/permutation_group.md) are described by the [Representation](documentation/symmetries/representation.md) object. As one-dimensional irreducible representations (irreps) are given by their characters $\chi(g)$ which is simply a list of (complex) numbers for every symmetry element $g$. A representation can be created by handing the group and the list of characters: 

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:symmetries_3"
	```
=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:symmetries_3"
	```

Upon creation of a representation, XDiag checks whether the group axioms as well as the homomorphism property of the characters is fulfilled, i.e.,

$$ f * g = h  \Rightarrow \chi(f) \cdot \chi(g) = \chi(h).$$

Representations can then be used to create symmetry adapted blocks, which can then be diagonalized independently:

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:symmetries_4"
	```
=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:symmetries_4"
	```
	
Also, permutation groups and representations can be defined in a TOML file and read in. Please refer to the [read_permutation_group](documentation/io/read_permutation_group.md) and [read_representation](documentation/io/read_representation.md) documentation on the specific format required. 

Using symmetries, thus, can give enhanced physical insights to the system, but also reduces the computational costs significantly. A common analysis tool for understanding quantum many-body systems is the **Tower of states** analysis. For an introduction, see e.g. 

> Studying Continuous Symmetry Breaking using Energy Level Spectroscopy
<br>
> Alexander Wietek, Michael Schuler, Andreas M. Läuchli
<br>
> arXiv:1704.08622 [cond-mat.str-el]<br>
> DOI: [10.48550/arXiv.1704.08622](https://doi.org/10.48550/arXiv.1704.08622)

## Further features

Some further features of XDiag include:

* Besides the [Spinhalf](documentation/blocks/spinhalf.md), further blocks include the [tJ](documentation/blocks/tJ.md) and [Electron](documentation/blocks/electron.md)

* Distributed blocks like [SpinhalfDistributed](documentation/blocks/spinhalf_distributed.md) and [tJDistributed](documentation/blocks/tJ_distributed.md) allow for distributed memory parallelization in C++.

* Real and imaginary time evolutions can be performed using the functions [time_evolve](documentation/algorithms/time_evolve.md), [imaginary_time_evolve](documentation/algorithms/imaginary_time_evolve.md).

* Operators can be symmetrized using the [symmetrize](documentation/operators/symmetrize.md) function. 

* Product states and random states can be created using the functions [product_state](documentation/states/create_state.md#product_state) and [random_state](documentation/states/create_state.md#random_state).

* Timing of certain parts of the code in C++ can be conventiently performed using the [Timing](documentation/utilities/timing.md) utilities.
