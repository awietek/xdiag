---
title: apply
---

Implements matrix multiplication of a sparse matrix $A$ in compressed-sparse-row (CSR) format with a vector $x$ or a matrix $X$,

$$ y = Ax, \quad Y = AX. $$

**Sources**<br>
[apply.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/algebra/sparse/apply.hpp)<br>
[apply.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/algebra/sparse/apply.cpp)<br>
[apply.jl](https://github.com/awietek/XDiag.jl/blob/main/src/algebra/sparse/apply.jl)

--- 

There are two interfaces to perform this operation.

1. The new vector $y$ or or matrix $Y$ is allocated internally and returned.

	=== "C++"
		```c++
		template <typename idx_t, typename coeff_t>
		arma::Col<coeff_t> apply(CSRMatrix<idx_t, coeff_t> const &A, arma::Col<coeff_t> const &x);
		template <typename idx_t, typename coeff_t>
		arma::Col<coeff_t> apply(CSRMatrix<idx_t, coeff_t> const &A, arma::Col<coeff_t> const &X);
		```
	=== "Julia"
		```julia
		apply(A::CSRMatrix{IdxT, CoeffT}, x::Vector{CoeffT})
		apply(A::CSRMatrix{IdxT, CoeffT}, X::Matrix{CoeffT})
		```

2. The return vector $y$ or matrix $Y$ has already been allocated, is given as an agrument to the apply function, where it is overwritten.

	=== "C++"
		```c++
		template <typename idx_t, typename coeff_t>
		void apply(CSRMatrix<idx_t, coeff_t> const &A,
                     arma::Col<coeff_t> const &x,
                     arma::Col<coeff_t> &y);
		template <typename idx_t, typename coeff_t>
		void apply(CSRMatrix<idx_t, coeff_t> const &A,
                     arma::Mat<coeff_t> const &X,
                     arma::Mat<coeff_t> &Y);
		```
	=== "Julia"
		```julia
		apply(A::CSRMatrix{IdxT, CoeffT}, x::Vector{CoeffT}, y::Vector{CoeffT})
		apply(A::CSRMatrix{IdxT, CoeffT}, X::Matrix{CoeffT}, Y::Matrix{CoeffT})
		```
---

## Parameters

| Name | Description                                                                                           |
|:-----|:------------------------------------------------------------------------------------------------------|
| A    | sparse matrix in compressed-sparse-row (CSR)format, see [sparse matrix types](sparse_matrix_types.md) |
| x, X | Input vector or matrix                                                                                |
| y, Y | Output vector or matrix                                                                               |

**Comment** XDiag only provides sparse matrix-vector and matrix-matrix multiplication functionality for CSR matrices, as these can be efficiently parallelized. When linking to the Intel MKL library in the C++ version, specialized sparse BLAS routines are called to improve performance.

---

## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:sparse_apply"
	```

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:sparse_apply"
	```
