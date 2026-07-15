---
title: Quick start
---

# Quick start

Welcome to XDiag! This quick start walks you through your very first exact diagonalization: computing the ground-state energy of the spin $S=1/2$ Heisenberg chain,
$$ H = J\sum_{\langle i,j \rangle} \mathbf{S}_i \cdot \mathbf{S}_j, $$
on a periodic one-dimensional lattice. Here $\mathbf{S}_i = (S_i^x, S_i^y, S_i^z)$ are the spin $S=1/2$ operators and $\langle i,j \rangle$ denotes summation over nearest-neighbor sites $i$ and $j$.

XDiag comes in two flavors that share the same API — a [Julia](#installation) package for interactive use and a [C++](#installation) library for maximum performance. Pick whichever you prefer below. For a thorough, step-by-step introduction, see the [User Guide](user_guide/index.md).

## Installation

=== "Julia"
	Enter the package mode in the Julia REPL by typing `]` and run
	```julia
	add XDiag
	```
	That's it — no compilation step is required.

=== "C++"
	First compile and install the XDiag library, following the [library compilation](documentation/compilation/compilation.md#library-compilation) instructions. The application code below is then compiled against it.

## Your first calculation

The complete program to set up the Heisenberg chain and compute its ground-state energy reads:

=== "Julia"
	```julia
	--8<-- "examples/spinhalf_chain_e0/main.jl"
	```
=== "C++"
	```c++
	--8<-- "examples/spinhalf_chain_e0/main.cpp"
	```

In just a few lines we define the Hilbert space (a `Spinhalf` block), build the Hamiltonian as an `OpSum`, and obtain the ground-state energy with `eigval0`.

!!! tip "C++ error traces"

	The `try / catch` clause implements an error-trace mechanism, activated by `error_trace`. While optional, we recommend it for every XDiag application, as it produces a readable traceback when a runtime error occurs.

!!! info "C++ compilation"

	Once written, the application is compiled with CMake — see the [application compilation](documentation/compilation/compilation.md#application-compilation) instructions. The C++ version also lets you optimize the compiled code for your target architecture, which can give a substantial speed-up (see [optimization](documentation/compilation/compilation.md#optimization)).

## Next steps

- Work through the [User Guide](user_guide/index.md) for a guided tour of all core features.
- Browse the [Example collection](examples.md) for fully worked applications.
- Consult the [Documentation](documentation/index.md) for the complete API reference.
