---
title: eigval0
---

Computes the groud state energy of a Hermitian operator on a block by using an iterative Lanczos algorithm. This function is a shortcut for the [eigvals_lanczos](eigvals_lanczos.md) function. We refer to [eigvals_lanczos](eigvals_lanczos.md) for further details on the algorithm and the convergence criterion.

The algorithm can be run either *on-the-fly* (matrix-free) or using a *sparse matrix* in the compressed-sparse-row format (see [CSRMatrix](../algebra/sparse/sparse_matrix_types.md)).

**Sources**<br>
[sparse_diag.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/algorithms/sparse_diag.hpp)<br>
[sparse_diag.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/algorithms/sparse_diag.cpp)<br>
[sparse_diag.jl](https://github.com/awietek/XDiag.jl/blob/main/src/algorithms/sparse_diag.jl)

---

## Definition

#### On-the-fly

=== "C++"
    ```c++
    double eigval0(OpSum const &ops, Block const &block, double precision = 1e-12,
                   int64_t max_iterations = 1000, int64_t random_seed = 42);
	```

=== "Julia"
	```julia
    eigval0(ops::OpSum, block::Block;	precision::Float64 = 1e-12, 
	        max_iterations::Int64 = 1000, random_seed::Int64 = 42)::Float64
	```
	
#### Sparse matrix

=== "C++"
    ```c++
    double eigval0(CSRMatrix<idx_t, coeff_t> const &ops,
		Block const &block, double precision = 1e-12,
		int64_t max_iterations = 1000, int64_t random_seed = 42);
	```

=== "Julia"
	```julia
    eigval0(ops::CSRMatrix, block::Block; precision::Float64 = 1e-12, 
		max_iterations::Int64 = 1000, random_seed::Int64 = 42)::Float64
	```


---

## Parameters

| Name           | Description                                                                                                           | Default |
|:---------------|:----------------------------------------------------------------------------------------------------------------------|---------|
| ops            | [OpSum](../operators/opsum.md) or [CSRMatrix](../algebra/sparse/sparse_matrix_types.md) defining a Hermitian operator |         |
| block          | block on which the operator is defined                                                                                |         |
| precision      | accuracy of the computed ground state                                                                                 | 1e-12   |
| max_iterations | maximum number of iterations                                                                                          | 1000    |
| random_seed    | random seed for setting up the initial vector                                                                         | 42      |

---

## Returns

| Type        | Description                      |
|:------------|:---------------------------------|
| real number | lowest lying eigenvalue of `ops` |

---

## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:eigval0"
	```
	
=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:eigval0"
	```


