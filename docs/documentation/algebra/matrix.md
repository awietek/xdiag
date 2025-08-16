---
title: matrix
---

Creates a numerical matrix with real (`matrix`) or complex (`matrixC`) coefficients given an [Op](../operators/op.md) or [OpSum](../operators/opsum.md) on a certain block. 

**Sources**<br>
[matrix.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/algebra/matrix.hpp)<br>
[matrix.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/algebra/matrix.cpp)<br>
[matrix.jl](https://github.com/awietek/XDiag.jl/blob/main/src/algebra/matrix.jl)

---

## Definition

A matrix can be created in two ways:

1. Only the input block on which the operator is defined is given. The output block is calculated and eventually created automatically.

	=== "C++"
		```c++
		arma::mat matrix(Op const &op, Block const &block);
		arma::mat matrix(OpSum const &ops, Block const &block);
		arma::cx_mat matrixC(Op const &op, Block const &block);
		arma::cx_mat matrixC(OpSum const &ops, Block const &block);
		```
	=== "Julia"
		```julia
		matrix(op::Op, block::Block)
		matrix(ops::OpSum, block::Block)
		```

**Comment:** In Julia, depending on whether a real/complex matrix is generated also a real/complex matrix is returned. The C++ version has to return a fixed type. If a real matrix is desired, use the functions `coo_matrix` or `coo_matrix_32`. If a complex matrix is desired, use the functions `coo_matrixC` and `coo_matrixC_32`.

---

## Parameters

| Name             | Description                                                                      |
|:-----------------|:---------------------------------------------------------------------------------|
| ops              | [OpSum](../operators/opsum.md) or [Op](../operators/op.md) defining the operator |
| block / block_in | input block on which the operator is defined                                     |
| block_out        | output block the operator maps the input block to                                |

---

## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:matrix"
	```

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:matrix"
	```
