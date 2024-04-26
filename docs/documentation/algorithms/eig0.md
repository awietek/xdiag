---
title: eig0
---

# eig0

```julia
eig0(bondlist, block; precision, max_iterations, force_complex, random_seed)
```

Computes the lowest lying eigenvalue of an operator.

##Parameters

| Name           | Description                                                            | Default |
|:---------------|:-----------------------------------------------------------------------|---------|
| bondlist       | BondList defining the bonds of the operator                            |         |
| block          | block on which the operator is defined                                 |         |
| precision      | accuracy of the computed ground state                                  | 1e-12   |
| max_iterations | maximum number of iterations                                           | 1000    |
| force_complex  | whether or not computation should be forced to have complex arithmetic | false   |
| random_seed    | random seed for setting up the initial vector                          | 42      |

##Returns

| Type        | Description             |
|:------------|:------------------------|
| real number | lowest lying eigenvalue |
| State       | groundstate             |

##Definition

=== "C++"

	```c++
	std::tuple<double, State>
	eig0(BondList const &bondlist, block_variant_t const &block,
		 double precision = 1e-12, int64_t max_iterations = 1000,
		 bool force_complex = false, int64_t random_seed = 42);
	```

=== "Julia"

	``` julia
	function eig0(bonds::BondList, block::Spinhalf;
	 	          precision::Real=1e-12, maxiter::Integer=1000,
                  force_complex::Bool=false, seed::Integer=42)
	```

Source: [sparse_diag.hpp](https://github.com/awietek/xdiag/blob/master/xdiag/algorithms/sparse_diag.hpp)
