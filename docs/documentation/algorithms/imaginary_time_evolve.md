---
title: imaginary_time_evolve
---

Computes the imaginary-time evolution, 

$$\vert \psi(\tau) \rangle = e^{-(H - \delta) \tau} \vert \psi_0\rangle,$$ 

of a [State](../states/state.md) $\vert \psi_0 \rangle$ and a Hermitian operator $H$ using an iterative algorithm. $\delta$ here denotes a real number which can be chosen as the ground state energy $\delta=E_0$ of $H$.

The algorithm can be run either *on-the-fly* (matrix-free) or using a *sparse matrix* in the compressed-sparse-row format (see [CSRMatrix](../algebra/sparse/sparse_matrix_types.md)).

**Sources** <br>
[imaginary_time_evolve.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/algorithms/imaginary_time_evolve.hpp)<br>
[imaginary_time_evolve.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/algorithms/imaginary_time_evolve.cpp)<br>
[imaginary_time_evolve.jl](https://github.com/awietek/XDiag.jl/blob/main/src/algorithms/imaginary_time_evolve.jl)
---

## Definition

#### On-the-fly

The method is provided in two variants:

1. Returning a new state while the input state remains untouched. This variant is safe to use and simple to code.

	=== "C++"

		```c++
		State imaginary_time_evolve(OpSum const &H, State psi0, double time,
                                    double precision = 1e-12, double shift = 0.);
		```
	=== "Julia"
		``` julia
		imaginary_time_evolve(ops::OpSum, psi0::State, time::Float64; precision::Float64 = 1e-12, 
		                      shift::Float64=0.0)::State
		```
		
2. An *inplace* variant `imaginary_time_evolve_inplace`, where the input state is overwritten and contains the time evolved state upon exit. This version is more memory efficient than `imaginary_time_evolve`.

	=== "C++"

		```c++
		void imaginary_time_evolve_inplace(OpSum const &H, State &psi0, double time,
                                           double precision = 1e-12, shift = 0.);
		```
	=== "Julia"
		``` julia
		imaginary_time_evolve_inplace(ops::OpSum, psi0::State, time::Float64; 
		                              precision::Float64 = 1e-12, shift::Float64=0.0)
		```
		
#### Sparse matrix

1. Returning a new state while the input state remains untouched. This variant is safe to use and simple to code.

	=== "C++"
		```c++
		template <typename idx_t, typename coeff_t>
		State imaginary_time_evolve(CSRMatrix<idx_t, coeff_t> const &H,
			State psi, double time, double precision = 1e-12, double shift = 0.);
		```
	=== "Julia"
		``` julia
		imaginary_time_evolve(H::CSRMatrix, psi0::State, time::Float64; 
			precision::Float64=1e-12, shift::Float64=0.0)::State
		```
		
2. An *inplace* variant `imaginary_time_evolve_inplace`, where the input state is overwritten and contains the time evolved state upon exit. This version is more memory efficient than `imaginary_time_evolve`.

	=== "C++"

		```c++
		template <typename idx_t, typename coeff_t>
		void imaginary_time_evolve_inplace(CSRMatrix<idx_t, coeff_t> const &H,
			State &psi, double time, double precision = 1e-12, double shift = 0.);
		```
	=== "Julia"
		``` julia
		imaginary_time_evolve_inplace(H::CSRMatrix, psi0::State, time::Float64; 
			precision::Float64 = 1e-12, shift::Float64=0.0)
		```

---

## Parameters

| Name      | Description                                                                                                                                    | Default |
|:----------|:-----------------------------------------------------------------------------------------------------------------------------------------------|---------|
| H         | [OpSum](../operators/opsum.md) or [CSRMatrix](../algebra/sparse/sparse_matrix_types.md) defining the hermitian operator $H$ for time evolution |         |
| psi0      | initial [State](../states/state.md) $\vert \psi_0 \rangle$ of the time evolution                                                               |         |
| time      | time $\tau$ until which the state is evolved                                                                                                   |         |
| precision | accuracy of the computed time evolved state $\vert \psi(\tau) \rangle$                                                                         | 1e-12   |
| shift     | the offset $\delta$ when computing $\vert \psi(t) \rangle = e^{-(H - \delta) \tau} \vert \psi_0\rangle$                                        | 0.0     |

The routine calls the subroutine [evolve_lanczos](evolve_lanczos.md) implementing a Lanczos algorithm to perform the evolution. This routine can also be called explicitly if more control is desired. Please also confer to the page [evolve_lanczos](evolve_lanczos.md) for further details on the specifics of the algorithm. The parameter `shift` can be used to turn all eigenvalues of the matrix $H - \delta \;\textrm{Id}$ positive whenever $\delta < E_0$, where $E_0$ denotes the ground state energy of $H$.

---

## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:imaginary_time_evolve"
	```
=== "Julia"
	```julia
	--8<-- "examples/usage_examples/main.jl:imaginary_time_evolve"
	```	
