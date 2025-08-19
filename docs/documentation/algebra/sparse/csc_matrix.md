---
title: csc_matrix
---

Creates a sparse matrix in the compressed-sparse-column (CSC) format. The sparse matrix can be constructed with real `csc_matrix` or complex `csc_matrixC` coefficients given an [OpSum](../../operators/opsum.md) on a certain block. 

By default the numbers describing the column and row indices are 64-bit signed integers. However, to save memory also 32-bit signed integer coded matrices can be created using the `csc_matrix_32` or `csc_matrixC_32` functions. Notice, that in this case the largest allowed row and column dimension is $2^{31} - 1= 2.147.483.647$.

For a description of the CSC sparse matrix format, see [Sparse matrix types](sparse_matrix_types.md).

**Sources**<br>
[csc_matrix.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/algebra/sparse/csc_matrix.hpp)<br>
[csc_matrix.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/algebra/sparse/csc_matrix.cpp)<br>
[csc_matrix.jl](https://github.com/awietek/XDiag.jl/blob/main/src/algebra/sparse/csc_matrix.jl)

---

## Definition

The sparse matrix can be created in two ways:

1. Only the input block on which the operator is defined is given. The output block is calculated and eventually created automatically.

	=== "C++"
		```c++
		CSCMatrix<int64_t, double> csc_matrix(OpSum const &ops, Block const &block, idx_t i0 = 0);
		CSCMatrix<int64_t, complex> csc_matrixC(OpSum const &ops, Block const &block, idx_t i0 = 0);
		CSCMatrix<int32_t, double> csc_matrix_32(OpSum const &ops, Block const &block, idx_t i0 = 0);
		CSCMatrix<int32_t, complex> csc_matrixC_32(OpSum const &ops, Block const &block, idx_t i0 = 0);
		```
	=== "Julia"
		```julia
		csc_matrix(ops::OpSum, block::Block, i0::Int64=1)
		csc_matrix_32(ops::OpSum, block::Block, i0::Int64=1)
		```
		
2. The output block is also handed as an argument. The compatibility of quantum numbers is checked. This way the output block is not created automatically and, thus, can be used to save computation time if the output block appears repeatedly in the computation.

	=== "C++"
		```c++
		CSCMatrix<int64_t, double> csc_matrix(OpSum const &ops, Block const &block_in, Block const &block_out, idx_t i0 = 0);
		CSCMatrix<int64_t, complex> csc_matrixC(OpSum const &ops, Block const &block_in, Block const &block_out, idx_t i0 = 0);
		CSCMatrix<int32_t, double> csc_matrix_32(OpSum const &ops, Block const &block_in, Block const &block_out, idx_t i0 = 0);
		CSCMatrix<int32_t, complex> csc_matrixC_32(OpSum const &ops, Block const &block_in, Block const &block_out, idx_t i0 = 0);
		```
	=== "Julia"
		```julia
		csc_matrix(ops::OpSum, block_in::Block, block_out::Block, i0::Int64=1)
		csc_matrix_32(ops::OpSum, block_in::Block, block_out::Block, i0::Int64=1)
		```
		
**Comment:** In Julia, depending on whether a real/complex matrix is generated also a real/complex matrix is returned. The C++ version has to return a fixed type. If a real matrix is desired, use the functions `csc_matrix` or `csc_matrix_32`. If a complex matrix is desired, use the functions `csc_matrixC` and `csc_matrixC_32`.

---

## Parameters

| Name             | Description                                                                        | Default            |
|:-----------------|:-----------------------------------------------------------------------------------|--------------------|
| ops              | [OpSum](../../operators/opsum.md) defining the operator                            |                    |
| block / block_in | input block on which the operator is defined                                       |                    |
| block_out        | output block the operator maps the input block to                                  |                    |
| i0               | integer saying whether integers are counted from 0 or 1, needs to be either 0 or 1 | 0 (C++), 1 (Julia) |

---

## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:csc_matrix"
	```

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:csc_matrix"
	```
