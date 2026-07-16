---
title: eigvals
---

Computes the `neigs` algebraically smallest eigenvalues of a hermitian operator using the [LOBPCG](eigs_lobpcg.md) block eigensolver. In contrast to the [Lanczos](eigvals_lanczos.md) algorithm, LOBPCG reliably resolves excited states and their degeneracies, at the cost of a higher memory footprint. Only the eigenvalues are returned; use [eigs](eigs.md) to additionally obtain the eigenvectors.

The algorithm can be run either *on-the-fly* (matrix-free) or using a *sparse matrix* in the compressed-sparse-row format (see [CSRMatrix](../kernels/sparse/sparse_matrix_types.md)).

**Sources:** [sparse_diag.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/linalg/sparse_diag.hpp) · [sparse_diag.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/linalg/sparse_diag.cpp)

## Definition

#### On-the-fly

=== "Julia"
	```julia
	eigvals(ops::OpSum, block::Block, neigs::Int64; precision::Float64 = 1e-12,
	        max_iterations::Int64 = 1000, random_seed::Int64 = 42)
	```
=== "C++"
	```c++
	arma::vec eigvals(OpSum const &ops, Block const &block, int64_t neigs,
	                  double precision = 1e-12, int64_t max_iterations = 1000,
	                  int64_t random_seed = 42);
	```

#### Sparse matrix

=== "Julia"
	```julia
	eigvals(ops::CSRMatrix, block::Block, neigs::Int64; precision::Float64 = 1e-12,
	        max_iterations::Int64 = 1000, random_seed::Int64 = 42)
	```
=== "C++"
	```c++
	template <typename idx_t, typename coeff_t>
	arma::vec eigvals(CSRMatrix<idx_t, coeff_t> const &ops, Block const &block,
	                  int64_t neigs, double precision = 1e-12,
	                  int64_t max_iterations = 1000, int64_t random_seed = 42);
	```

## Parameters

| Name           | Description                                                                                                                | Default |
|:---------------|:---------------------------------------------------------------------------------------------------------------------------|---------|
| ops            | [OpSum](../operators/opsum.md) or [CSRMatrix](../kernels/sparse/sparse_matrix_types.md) defining the bonds of the operator |         |
| block          | block on which the operator is defined                                                                                     |         |
| neigs          | number of (lowest) eigenvalues to compute                                                                                  |         |
| precision      | accuracy of the computed eigenvalues                                                                                       | 1e-12   |
| max_iterations | maximum number of iterations                                                                                               | 1000    |
| random_seed    | random seed for setting up the initial block of vectors                                                                    | 42      |

## Returns

An `arma::vec` (Julia: `Vector{Float64}`) holding the `neigs` lowest eigenvalues in ascending order.

## Usage Example

=== "Julia"
	```julia
	block = Spinhalf(12)
	ops = OpSum()
	for i in 1:12
	    ops += "J" * Op("SdotS", [i, mod1(i + 1, 12)])
	end
	ops["J"] = 1.0
	eigenvalues = eigvals(ops, block, 3)
	```
=== "C++"
	```c++
	auto block = Spinhalf(12);
	auto ops = OpSum();
	for (int i = 0; i < 12; ++i) {
	  ops += "J" * Op("SdotS", {i, (i + 1) % 12});
	}
	ops["J"] = 1.0;
	arma::vec eigenvalues = eigvals(ops, block, 3);
	```
