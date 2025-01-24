---
title: eig0
---

Computes the groud state energy and the ground state of a Hermitian operator on a block by using an iterative Lanczos algorithm.

**Sources** [sparse_diag.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/algorithms/sparse_diag.hpp) [sparse_diag.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/algorithms/sparse_diag.cpp)

## Definition

=== "C++"

    ```c++
	std::tuple<double, State> eig0(OpSum const &ops, Block const &block,
		double precision = 1e-12,
		int64_t max_iterations = 1000,
		bool force_complex = false,
		int64_t random_seed = 42);
	```

=== "Julia"
	
	``` julia
    function eig0(
		ops::OpSum,
		block::Block;
		precision::Real = 1e-12,
		maxiter::Int64 = 1000,
		force_complex::Bool = false,
		seed::Int64 = 42,
	)
	```



## Parameters

| Name           | Description                                                            | Default |
|:---------------|:-----------------------------------------------------------------------|---------|
| ops            | [OpSum](../operators/opsum.md) defining the bonds of the operator      |         |
| block          | block on which the operator is defined                                 |         |
| precision      | accuracy of the computed ground state                                  | 1e-12   |
| max_iterations | maximum number of iterations                                           | 1000    |
| force_complex  | whether or not computation should be forced to have complex arithmetic | false   |
| random_seed    | random seed for setting up the initial vector                          | 42      |

## Returns

| Type        | Description             |
|:------------|:------------------------|
| real number | lowest lying eigenvalue |
| State       | groundstate             |

## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:eig0"
	```
	
=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:eig0"
	```

