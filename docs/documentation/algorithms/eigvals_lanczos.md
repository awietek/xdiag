---
title: eigvals_lanczos
---

Performs an iterative eigenvalue calculation using the Lanczos algorithm. Returns the tridiagonal matrix, eigenvalues, number of iterations and the stopping criterion.

The algorithm can be run either *on-the-fly* (matrix-free) or using a *sparse matrix* in the compressed-sparse-row format (see [CSRMatrix](../algebra/sparse/sparse_matrix_types.md)).

**Sources**<br>
[eigvals_lanczos.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/algorithms/lanczos/eigvals_lanczos.hpp)<br>
[eigvals_lanczos.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/algorithms/lanczos/eigvals_lanczos.cpp)<br>
[eigvals_lanczos.jl](https://github.com/awietek/XDiag.jl/blob/main/src/algorithms/lanczos/eigvals_lanczos.jl)

---

## Definition

#### On-the-fly

The Lanczos algorithm can be run in thre distinct ways:

1. A random intial state $|\psi_0\rangle = |r\rangle$ with normal distributed entries is used.

	=== "C++"

		```c++
		EigvalsLanczosResult
		eigvals_lanczos(OpSum const &ops, Block const &block, int64_t neigvals = 1,
                    	double precision = 1e-12, int64_t max_iterations = 1000,
                        double deflation_tol = 1e-7, int64_t random_seed = 42);
		```
	
	=== "Julia"

		``` julia
		eigvals_lanczos(ops::OpSum, block::Block; neigvals::Int64 = 1, 
		                precision::Float64 = 1e-12,	max_iterations::Int64 = 1000,
                        deflation_tol::Float64 = 1e-7, random_seed::Int64 = 42)
		```

2. The initial state $|\psi_0\rangle$ is explicitly specified. 

	=== "C++"
		```c++
		EigvalsLanczosResult 
		eigvals_lanczos(OpSum const &ops, State psi0, int64_t neigvals = 1,
	                    double precision = 1e-12, int64_t max_iterations = 1000,
						double deflation_tol = 1e-7);
     	```
	=== "Julia"

		``` julia
		eigvals_lanczos(ops::OpSum, psi0::State; neigvals::Int64 = 1, 
		                precision::Float64 = 1e-12,	max_iterations::Int64 = 1000,
                        deflation_tol::Float64 = 1e-7)
		```
		
	Notice this version copies the initial state, which requires memory but keeps the orginal state intact.

3. The initial state $|\psi_0\rangle$ is explicitly specified and **overwritten** in the process. This version can save memory, but the initial state  $|\psi_0\rangle$ cannot be used later.

	=== "C++"
		```c++
		EigvalsLanczosResult 
		eigvals_lanczos_inplace(OpSum const &ops, State &psi0, int64_t neigvals = 1,
	                        	double precision = 1e-12, int64_t max_iterations = 1000,
                                double deflation_tol = 1e-7);
     	```
	=== "Julia"

		``` julia
		eigvals_lanczos_inplace(ops::OpSum, psi0::State; neigvals::Int64 = 1, 
		                        precision::Float64 = 1e-12,	max_iterations::Int64 = 1000,
                                deflation_tol::Float64 = 1e-7)
		```

#### Sparse matrix

1. A random intial state $|\psi_0\rangle = |r\rangle$ with normal distributed entries is used.

	=== "C++"

		```c++
		template <typename idx_t, typename coeff_t>
		EigvalsLanczosResult
		eigvals_lanczos(CSRMatrix<idx_t, coeff_t> const &ops, Block const &block,
                int64_t neigvals = 1, double precision = 1e-12,
                int64_t max_iterations = 1000, double deflation_tol = 1e-7,
                int64_t random_seed = 42);
		```
	
	=== "Julia"

		``` julia
		eigvals_lanczos(ops::CSRMatrix, block::Block; neigvals::Int64 = 1, 
			precision::Float64 = 1e-12, max_iterations::Int64 = 1000,
			deflation_tol::Float64 = 1e-7, random_seed::Int64 = 42)::EigvalsLanczosResult
		```

2. The initial state $|\psi_0\rangle$ is explicitly specified. 

	=== "C++"
		```c++
		template <typename idx_t, typename coeff_t>
	    EigvalsLanczosResult
		eigvals_lanczos(CSRMatrix<idx_t, coeff_t> const &ops, State psi0,
                int64_t neigvals = 1, double precision = 1e-12,
                int64_t max_iterations = 1000, double deflation_tol = 1e-7);
     	```
	=== "Julia"

		``` julia
	    eigvals_lanczos(ops::CSRMatrix, psi0::State; neigvals::Int64 = 1, 
			precision::Float64 = 1e-12, max_iterations::Int64 = 1000,
			deflation_tol::Float64 = 1e-7)::EigvalsLanczosResult
		```
		
	Notice this version copies the initial state, which requires memory but keeps the orginal state intact.

3. The initial state $|\psi_0\rangle$ is explicitly specified and **overwritten** in the process. This version can save memory, but the initial state  $|\psi_0\rangle$ cannot be used later.

	=== "C++"
		```c++
		template <typename idx_t, typename coeff_t>
		EigvalsLanczosResult 
		eigvals_lanczos_inplace(CSRMatrix<idx_t, coeff_t> const &ops, 
			State &psi0, int64_t neigvals = 1, double precision = 1e-12,
			int64_t max_iterations = 1000, double deflation_tol = 1e-7);
     	```
	=== "Julia"

		``` julia
		eigvals_lanczos_inplace(ops::CSRMatrix, psi0::State; neigvals::Int64 = 1, 
			precision::Float64 = 1e-12, max_iterations::Int64 = 1000,
			deflation_tol::Float64 = 1e-7)::EigvalsLanczosResult
		```


---

## Parameters

| Name           | Description                                                                                                                | Default |
|:---------------|:---------------------------------------------------------------------------------------------------------------------------|---------|
| ops            | [OpSum](../operators/opsum.md) or [CSRMatrix](../algebra/sparse/sparse_matrix_types.md) defining the bonds of the operator |         |
| block          | block on which the operator is defined                                                                                     |         |
| psi0           | Initial [State](../states/state.md) from which the Lanczos algorithm is started                                            |         |
| neigvals       | number $k$ of eigenvalue to converge                                                                                       | 1       |
| precision      | accuracy of the computed ground state                                                                                      | 1e-12   |
| max_iterations | maximum number of iterations                                                                                               | 1000    |
| deflation_tol  | tolerance for deflation, i.e. breakdown of Lanczos due to Krylow space exhaustion                                          | 1e-7    |
| random_seed    | random seed for setting up the initial vector                                                                              | 42      |

---

## Returns

A struct of type `EigvalsLanczosResult` with the following entries.

| Entry       | Description                                             |
|:------------|:--------------------------------------------------------|
| alphas      | diagonal elements of the tridiagonal matrix             |
| betas       | off-diagonal elements of the tridiagonal matrix         |
| eigenvalues | the computed Ritz eigenvalues of the tridiagonal matrix |
| niterations | number of iterations performed                          |
| criterion   | string denoting the reason why the algorithm stopped    |

---

## Convergence criterion

The algorithm terminates if the $k$-th ($k$ is the argument `neigvals`) approximate eigenvalue changes only by a fraction smaller than $\epsilon$ ($k$ is the argument `precision`), i.e.

$$ (\tilde{e}_k^{(n)} - \tilde{e}_k^{(n-1)}) / \tilde{e}_k^{(n)} < \epsilon.$$

Here, $\tilde{e}_k^{(n)}$ denotes the Lanczos approximation to the $k$-th eigenvalue after $n$ iterations.

---

## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:eigvals_lanczos"
	```

=== "Julia"
	```julia
	--8<-- "examples/usage_examples/main.jl:eigvals_lanczos"
	```
