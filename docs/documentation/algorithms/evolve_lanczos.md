---
title: evolve_lanczos
---

Computes the exponential of a Hermitian operator $H$ with an arbitrary real or complex prefactor $z$ applied to a [State](../states/state.md) $\vert \psi_0\rangle$, 

$$\vert \psi(z) \rangle = e^{z(H - \delta)} \vert \psi_0\rangle.$$ 

Here, $\delta$ denotes a real number shifting the spectrum of $H$. The algorithm implemented is described in the following publication.

> On Krylov Subspace Approximations to the Matrix Exponential Operator<br>
> Marlis Hochbruck and Christian Lubich<br>
> SIAM Journal on Numerical Analysis, Vol. 34, Iss. 5 (1997)<br>
> DOI: [10.1137/S0036142995280572](https://doi.org/10.1137/S0036142995280572)

**Sources**<br>
[evolve_lanczos.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/algorithms/evolve_lanczos.hpp)<br>
[evolve_lanczos.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/algorithms/evolve_lanczos.cpp)<br>
[evolve_lanczos.jl](https://github.com/awietek/XDiag.jl/blob/main/src/algorithms/evolve_lanczos.jl)

---

## Definition

The method is provided in two variants:

1. Returning a new state while the input state remains untouched. This variant is safe to use and simple to code.

	=== "C++"
		```c++
		EvolveLanczosResult
		evolve_lanczos(OpSum const &H, State psi, double t, double precision = 1e-12,
      		           double shift = 0., bool normalize = false,
                       int64_t max_iterations = 1000, double deflation_tol = 1e-7);

		EvolveLanczosResult
		evolve_lanczos(OpSum const &H, State psi, complex z, double precision = 1e-12,
      		           double shift = 0., bool normalize = false,
                       int64_t max_iterations = 1000, double deflation_tol = 1e-7);
		```
		
	=== "Julia"
		```julia
		evolve_lanczos(H::OpSum, psi::State, t::Float64; precision::Float64 = 1e-12,
      		           shift::Float64=0.0, normalize::Bool=false,
                       max_iterations::Int64 = 1000, deflation_tol::Float64 = 1e-7)

		evolve_lanczos(H::OpSum, psi::State, z::ComplexF64; precision::Float64 = 1e-12,
	                   shift::Float64=0.0, normalize::Bool=false,
				       max_iterations::Int64 = 1000, deflation_tol::Float64 = 1e-7)
		```
		

2. An *inplace* variant `evolve_lanczos_inplace`, where the input state is overwritten and contains the time evolved state upon exit. This version is more memory efficient than `evolve_lanczos`.

	=== "C++"
		```c++
		EvolveLanczosInplaceResult
		evolve_lanczos_inplace(OpSum const &H, State &psi, double t, 
		                       double precision = 1e-12, double shift = 0.,
							   bool normalize = false, int64_t max_iterations = 1000, 
							   double deflation_tol = 1e-7);

		EvolveLanczosInplaceResult
		evolve_lanczos_inplace(OpSum const &H, State &psi, complex z, 
		                       double precision = 1e-12, double shift = 0.,
							   bool normalize = false, int64_t max_iterations = 1000, 
							   double deflation_tol = 1e-7);
		```
	=== "Julia"
		```julia
		evolve_lanczos_inplace(H::OpSum, psi::State, t::Float64; precision::Float64 = 1e-12,
   	                           shift::Float64=0.0, normalize::Bool=false,
                               max_iterations::Int64 = 1000, deflation_tol::Float64 = 1e-7)

		evolve_lanczos_inplace(H::OpSum, psi::State, z::ComplexF64; precision::Float64 = 1e-12,
	                           shift::Float64=0.0, normalize::Bool=false,
				               max_iterations::Int64 = 1000, deflation_tol::Float64 = 1e-7)
		```

---

## Parameters

| Name           | Description                                                                                             | Default |
|:---------------|:--------------------------------------------------------------------------------------------------------|---------|
| H              | [OpSum](../operators/opsum.md) defining the hermitian operator $H$ for time evolution                   |         |
| psi0           | initial [State](../states/state.md) $\vert \psi_0 \rangle$ of the time evolution                        |         |
| time           | time $\tau$ until which the state is evolved                                                            |         |
| precision      | accuracy of the computed time evolved state $\vert \psi(t) \rangle$                                     | 1e-12   |
| shift          | the offset $\delta$ when computing $\vert \psi(t) \rangle = e^{-(H - \delta) \tau} \vert \psi_0\rangle$ | 0.0     |
| normalize      | flag whether or not the evolved state should be normalized                                              | false   |
| max_iterations | maximum number of Lanczos iterations performed                                                          | 1000    |
| deflation_tol  | tolerance for deflation, i.e. breakdown of Lanczos due to Krylow space exhaustion                       | 1e-7    |

The parameter `shift` can be used to turn all eigenvalues of the matrix $H - \delta \;\textrm{Id}$ positive whenever $\delta < E_0$, where $E_0$ denotes the ground state energy of $H$.

---

## Returns

A struct with the following entries

| Entry       | Description                                                                                       |
|:------------|:--------------------------------------------------------------------------------------------------|
| alphas      | diagonal elements of the Lanczos tridiagonal matrix                                               |
| betas       | off-diagonal elements of the Lanczos tridiagonal matrix                                           |
| eigenvalues | the computed Ritz eigenvalues of the tridiagonal matrix                                           |
| niterations | number of iterations performed                                                                    |
| criterion   | string denoting the reason why the Lanczosalgorithm stopped                                       |
| state       | time-evolved [State](../states/state.md) $\vert \psi(t)\rangle$ (not defined for inplace variant) |

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
	--8<-- "examples/usage_examples/main.cpp:evolve_lanczos"
	```
=== "Julia"
	```julia
	--8<-- "examples/usage_examples/main.jl:evolve_lanczos"
	```
		
