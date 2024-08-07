---
title: matrix
---

```julia
matrix(ops, block; force_complex=false)
matrixC(ops, block)  # c++ only
```

Creates the full matrix representation of a given [OpSum](../operators/opsum.md) on a block.

In Julia, depending on whether a real/complex matrix is generated also a  real/complex matrix is returned. The C++ version has to return a fixed type. If a real matrix is desired, use the function `matrix`. If a complex matrix is desired, use the function `matrixC`.

**Source** [matrix.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/algebra/matrix.hpp)

## Parameters

| Name          | Description                                                               |   |
|:--------------|:--------------------------------------------------------------------------|---|
| ops           | [OpSum](../operators/opsum.md) defining the operator                      |   |
| block         | block on which the operator is defined                                    |   |
| force_complex | flag to determine if returned matrix is forced to be complex (Julia only) |   |

## Returns
		
| Type   | Description                                |   |
|:-------|:-------------------------------------------|---|
| matrix | matrix representation of opsum on block |   |


## Definition

=== "Julia"
	```julia
    matrix(ops::Opsum, block::Block; force_complex::Bool=false)
	```

=== "C++"
	```c++
	template <typename block_t>
    arma::mat matrix(OpSum const &ops, block_t const &block);

    template <typename block_t>
    arma::cx_mat matrixC(OpSum const &ops, block_t const &block);
	```

## Usage Example

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:matrix"
	```

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:matrix"
	```
	
