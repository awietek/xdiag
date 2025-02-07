---
title: time_evolve_expokit
---

Computes the real-time evolution, 

$$\vert \psi(t) \rangle = e^{-iHt} \vert \psi_0\rangle,$$ 

of a [State](../states/state.md) $\vert \psi_0 \rangle$ and a Hermitian operator $H$ using the iterative algorithm implemented by [Expokit](https://www.maths.uq.edu.au/expokit/)

> Expokit: A Software Package for Computing Matrix Exponentials<br>
> Roger B. Sidje<br>
> ACM Trans. Math. Softw., 24(1):130-156, 1998. (1998)<br>
> DOI: [10.1145/285861.285868](https://doi.org/10.1145/285861.285868)

The algorithm features automatic stepsize control and computes approximate solutions with high precision according to our tests. Yet, the [evolve_lanczos](evolve_lanczos.md) implementation is currently faster and more memory efficient. 

**Sources**<br> 
[time_evolve_expokit.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/algorithms/time_evolution/time_evolve_expokit.hpp)<br>
[time_evolve_expokit.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/algorithms/time_evolution/time_evolve_expokit.cpp)<br>
[time_evolve_expokit.jl](https://github.com/awietek/XDiag.jl/blob/main/src/algorithms/time_evolution/time_evolve_expokit.jl)

---

## Definition

The method is provided in two variants:

1. Returning a new state while the input state remains untouched. This variant is safe to use and simple to code.

	=== "C++"
		```c++
		time_evolve_expokit_return_t time_evolve_expokit(
			OpSum const &ops, State state, double time, double precision = 1e-12,
			int64_t m = 30, double anorm = 0., int64_t nnorm = 2);
		```
	=== "Julia"
		```julia
		time_evolve_expokit(ops::OpSum, state::State, time::Float64;
		                    precision::Float=1e-12, m::Int64 = 30, 
							anorm::Float64 = 0.0, nnorm::Int64 = 2)
		```

2. An *inplace* variant `time_evolve_expokit_inplace`, where the input state is overwritten and contains the time evolved state upon exit. This version is more memory efficient than `time_evolve_expokit`.

	=== "C++"
		```c++
		time_evolve_expokit_inplace_return_t time_evolve_expokit_inplace(
			OpSum const &ops, State &state, double time, double precision = 1e-12,
			int64_t m = 30, double anorm = 0., int64_t nnorm = 2);
		```
	=== "Julia"
		```julia
		time_evolve_expokit_inplace(ops::OpSum, state::State, time::Float64;
		                            precision::Float=1e-12, m::Int64 = 30, 
							        anorm::Float64 = 0.0, nnorm::Int64 = 2)
		```
---

## Parameters

| Name      | Description                                                                           | Default |
|:----------|:--------------------------------------------------------------------------------------|---------|
| H         | [OpSum](../operators/opsum.md) defining the hermitian operator $H$ for time evolution |         |
| psi0      | initial [State](../states/state.md) $\vert \psi_0 \rangle$ of the time evolution      |         |
| time      | time $t$ until which the state is evolved                                             |         |
| precision | accuracy of the computed time evolved state                                           | 1e-12   |
| m         | dimension of used Krylov space, main memory requirement                               | 30      |
| anorm     | 1-norm estimate of the operator $H$, if unknown default 0. computes it fresh          | 0.      |
| nnorm     | number of random samples to estimate 1-norm, usually not more than 2 required         | 2       |

---

## Returns

A struct with the following entries

| Entry | Description                                                                                       |
|:------|:--------------------------------------------------------------------------------------------------|
| error | the computed error estimate during evolution                                                      |
| hump  | the "hump" as defined in Expokit [10.1145/285861.285868](https://doi.org/10.1145/285861.285868)   |
| state | time-evolved [State](../states/state.md) $\vert \psi(t)\rangle$ (not defined for inplace variant) |

---

## Convergence criterion

The algorithm is estimating the following error,

$$ \varepsilon = \parallel \vert \tilde{\psi}(t)\rangle - e^{z(H - \delta)} \vert \psi_0\rangle \parallel_2, $$

where $\vert \tilde{\psi}(t) \rangle$ denotes the approximation computed during the algorithm. As the exact solution is not available this error is estimated using the method described by *Algorithm 2* in

> Expokit: A Software Package for Computing Matrix Exponentials<br>
> Roger B. Sidje<br>
> ACM Trans. Math. Softw., 24(1):130-156, 1998. (1998)<br>
> DOI: [10.1145/285861.285868](https://doi.org/10.1145/285861.285868)

---

## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:time_evolve_expokit"
	```
	
=== "Julia"
	```julia
	--8<-- "examples/usage_examples/main.jl:time_evolve_expokit"
	```
	
