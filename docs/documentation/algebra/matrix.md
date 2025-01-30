---
title: matrix
---

Creates a numerical matrix with real (`matrix`) or complex (`matrixC`) coefficients given an [Op](../operators/op.md) or [OpSum](../operators/opsum.md) on a certain block. 

**Sources** [matrix.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/algebra/matrix.hpp), [matrix.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/algebra/matrix.cpp)

---

## Definition

A matrix can be created in two ways:

1. Only the input block on which the operator is defined is given. The output block is calculated and eventually created automatically.

	=== "C++"
		```c++
		template <class block_t>
		arma::mat matrix(Op const &op, block_t const &block);

		template <class block_t>
		arma::mat matrix(OpSum const &ops, block_t const &block);

		template <class block_t>
		arma::cx_mat matrixC(Op const &op, block_t const &block);

		template <typename block_t>
		arma::cx_mat matrixC(OpSum const &ops, block_t const &block);
		```
	=== "Julia"
		```julia
		matrix(op::Op, block::Block)
		matrix(ops::OpSum, block::Block)
		```

2. The output block is also handed as an argument. The compatibility of quantum numbers is checked. This way the output block is not created automatically and, thus, can be used to save computation time if the output block appears repeatedly in the computation.

	=== "C++"
		```c++
		template <class block_t>
		arma::mat matrix(Op const &op, block_t const &block_in, 
		                 block_t const &block_out);

		template <class block_t>
		arma::mat matrix(OpSum const &ops, block_t const &block_in, 
		                 block_t const &block_out);

		template <class block_t>
		arma::cx_mat matrixC(Op const &op, block_t const &block_in, 
		                     block_t const &block_out);

		template <typename block_t>
		arma::cx_mat matrixC(OpSum const &ops, block_t const &block_in, 
		                     block_t const &block_out);
		```
	=== "Julia"
		```julia
		matrix(op::Op, block_in::Block, block_out::Block)
		matrix(ops::OpSum, block_in::Block, block_out::Block)
		```

**Comment:** In Julia, depending on whether a real/complex matrix is generated also a  real/complex matrix is returned. The C++ version has to return a fixed type. If a real matrix is desired, use the function `matrix`. If a complex matrix is desired, use the function `matrixC`.

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
