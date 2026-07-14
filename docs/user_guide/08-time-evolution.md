---
title: Time evolution
---

# Time evolution

Besides diagonalization, XDiag also offers iterative algorithms to perform real- and imaginary-time evolutions of the form,

$$ |\phi(t)\rangle = e^{-iHt} |\psi_0\rangle \quad \text{or} \quad |\eta(\tau)\rangle = e^{-\tau H} |\psi_0\rangle. $$

Like the iterative diagonalization routines, these algorithms are *matrix-free*: they take the Hamiltonian as an [OpSum](../documentation/operators/opsum.md) and act on a [State](../documentation/states/state.md), without ever building a dense matrix.

## Real-time evolution

A real-time evolution $e^{-iHt}|\psi_0\rangle$ is performed with the [time_evolve](../documentation/linalg/time_evolve.md) function.

=== "Julia"
	```julia
	--8<-- "examples/user_guide/main.jl:usage_guide_iter3"
	```
=== "C++"
	```c++
	--8<-- "examples/user_guide/main.cpp:usage_guide_iter3"
	```

Two distinct algorithms are available. The first is the memory-efficient Lanczos algorithm described by [Hochbruck and Lubich (1997)](https://epubs.siam.org/doi/10.1137/S0036142995280572), which runs a Lanczos iteration twice, first to compute the tridiagonal matrix and then to build the time-evolved state. The second is the algorithm proposed in [Expokit](https://dl.acm.org/doi/10.1145/285861.285868). While the Expokit algorithm is computationally more efficient and highly accurate, it has higher memory requirements. The algorithm is selected through the optional `algorithm` argument, which defaults to `"lanczos"`; setting it to `"expokit"` selects the Expokit algorithm instead.

## Imaginary-time evolution

An imaginary-time evolution $e^{-\tau H}|\psi_0\rangle$ is performed with the [imaginary_time_evolve](../documentation/linalg/imaginary_time_evolve.md) function.

=== "Julia"
	```julia
	--8<-- "examples/user_guide/main.jl:usage_guide_time_imag"
	```
=== "C++"
	```c++
	--8<-- "examples/user_guide/main.cpp:usage_guide_time_imag"
	```

Since the operator $e^{-\tau H}$ amplifies the components of low-lying eigenstates exponentially, an optional `shift` argument is provided. It shifts the Hamiltonian $H \to H - \text{shift}$ during the evolution, which does not change the normalized state but avoids numerical overflow. A natural choice is to set `shift` to (an estimate of) the ground state energy.

## Fine-grained control

More direct control over both algorithms is provided by the functions [evolve_lanczos](../documentation/linalg/evolve_lanczos.md), [time_evolve_expokit](../documentation/linalg/time_evolve_expokit.md), and their in-place variants. These allow additional parameters to be set and return further information about the run. [evolve_lanczos](../documentation/linalg/evolve_lanczos.md) returns the coefficients of the tridiagonal Lanczos matrix (and accepts a real or complex time argument, covering both real- and imaginary-time evolution), while [time_evolve_expokit](../documentation/linalg/time_evolve_expokit.md) returns error estimates of the computed evolution. For a precise description of these methods, we refer to their respective documentations.
