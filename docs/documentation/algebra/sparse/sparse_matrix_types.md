---
title: Sparse matrix types
---

As opposed to dense matrices, there exist several distinct formats in which sparse matrices are stored. XDiag offers capabilities to create several common sparse matrix formats which can be applied for the provided iterative algorithms. An overview of common sparse matrix formats formats can also be found on [Wikipedia](https://en.wikipedia.org/wiki/Sparse_matrix), where also a detailed description and examples can be found. In the following we describe the three supported sparse matrix formats and how these are represented in the code.

## Coordinate (COO) format

The coordinate format is a simple format which consists of two integer arrays `nrows` and `ncols` which store the row and column indices of the non-zero elements, and an additional array `data` storing the actual numerical matrix entries as a real or complex number. Within XDiag, we represent this format as a simple struct, without any class functionality. 

=== "C++"
	```c++
	template <typename idx_t, typename coeff_t> 
	struct COOMatrix {
		idx_t nrows;             // number of rows in the matrix
		idx_t ncols;             // number of columns in the matrix
		arma::Col<idx_t> row;    // row indices
		arma::Col<idx_t> col;    // column indices
		arma::Col<coeff_t> data; // data entries
		idx_t i0;                // zero index, either 0 or 1
		bool ishermitian;        // flag whether matrix is hermitian/symmetric
	};
	```
=== "Julia"
	```julia
	struct COOMatrix{IdxT<:Integer,CoeffT<:Number}
		nrows::IdxT
		ncols::IdxT
		row::Vector{IdxT}
		col::Vector{IdxT} 
		data::Vector{CoeffT}
		i0::IdxT
		ishermitian::Bool
	end
	```

Notice, that in addition to the three arrays `row`, `col` `data`, also the total number of rows `nrows` and columns `ncols` is stored. Also, the entry `i0` is used in XDiag to indicate whether the matrix is 0-indexed (`i0=0`) or 1-indexed (`i0=1`). Lastly, whenever a sparse matrix is created in XDiag the flag `ishermitian` is set, which indicates whether or not the matrix is Hermitian. 

The COO format is simple to understand and can be created very efficiently using XDiag. The entries in the `row` and `col` arrays are not ordered in any specific manner. This format is a common format in many sparse matrix libraries to create more advances storage types. A disavantage of this format is the memory overhead by storing both row and column indices, which can for example be reduces in the compressed-sparse-row (CSR) or compressed-sparse-column (CSC) formats. Also, typically implementations of matrix-vector multiplications are often not efficient for the COO format, due to a lack of locality in memory.

To create matrices in the COO format, the functions [coo_matrix](coo_matrix.md) and [coo_matrix_32](coo_matrix.md) can be used.

## Compressed-sparse-row (CSR) format

The CSR format compresses information about the indices, by only storing the column indices in an integer array `col` and a "pointer" array `rowptr` (simply an integer index) to the first element of each row. Finally, the coefficients of the matrix are then stored in an array called `data`. For a detailed explanation of the format we refer to [Wikipedia](https://en.wikipedia.org/wiki/Sparse_matrix). Within XDiag, CSR matrices are represented by simple structs without any additional functionality. 

=== "C++"
	```c++
	template <typename idx_t, typename coeff_t> 
	struct CSRMatrix {
		idx_t nrows;             // number of rows in the matrix
		idx_t ncols;             // number of columns in the matrix
		arma::Col<idx_t> rowptr; // pointer to elements of a row (size: nrows+1)
		arma::Col<idx_t> col;    // columns for each row consecutively
		arma::Col<coeff_t> data; // data entries each row consecutively
		idx_t i0;                // zero index, either 0 or 1
		bool ishermitian;        // flag whether matrix is hermitian/symmetric
	};
	```
=== "Julia"
	```julia
	struct CSRMatrix{IdxT<:Integer,CoeffT<:Number}
		nrows::IdxT
		ncols::IdxT
		rowptr::Vector{IdxT}
		col::Vector{IdxT} 
		data::Vector{CoeffT}
		i0::IdxT
		ishermitian::Bool
	end
	```

Again the entries `nrows` and `ncols` denote the number of rows and columns, `i0` denotes whether the matrix is 0-indexed or 1-indexed, and the flag `ishermitian` is `true` if the sparse matrix is Hermitian.

CSR matrices created by XDiag are sorted in the sense that subsequent column indices within each row are monotonically increasing. CSR matrices (as well as CSC matrices) use less memory than the corresponding matrix in the COO format, as the information about row indices is compressed. Importantly, CSR matrix-vector or matrix-matrix multiplications can be efficiently parallelized, which makes this format highly attractive. The internal iterative algorithms provide support only for the CSR format, since with COO and CSC format it is very difficult to achieve comparable parallel performance as with CSR matrices. 

To create matrices in the CSR format, the functions [csr_matrix](csr_matrix.md) and [csr_matrix_32](csr_matrix.md) can be used.

## Compressed-sparse-column (CSC) format

The CSC format is similar to the CSR format, but with the role of columns and rows switched. Hence, there is an integer array `row` storing al the row indices and the `colptr` which stores the index of the first element in each column. The `data` array then again stores values of the matrix entries. A CSC-format matrix is represented by simple structs in XDiag:

=== "C++"
	```c++
	template <typename idx_t, typename coeff_t> 
	struct CSCMatrix {
		idx_t nrows;             // number of rows in the matrix
		idx_t ncols;             // number of columns in the matrix
		arma::Col<idx_t> colptr; // pointer to elements of a row (size: ncols+1)
		arma::Col<idx_t> row;    // columns for each row consecutively
		arma::Col<coeff_t> data; // data entries each row consecutively
		idx_t i0;                // zero index, either 0 or 1
		bool ishermitian;        // flag whether matrix is hermitian/symmetric
	};
	```
=== "Julia"
	```julia
	struct CSCMatrix{IdxT<:Integer,CoeffT<:Number}
		nrows::IdxT
		ncols::IdxT
		colptr::Vector{IdxT} 
		row::Vector{IdxT}
		data::Vector{CoeffT}
		i0::IdxT
		ishermitian::Bool
	end
	```

Again the entries `nrows` and `ncols` denote the number of rows and columns, `i0` denotes whether the matrix is 0-indexed or 1-indexed, and the flag `ishermitian` is `true` if the sparse matrix is Hermitian.

CSC matrices are the standard choice for storing sparse matrices, both in Julia and armadillo. While this format offers several advantages like fast-column access, parallelizations of matrix-vector multiplications are more difficult, see for example this [discussion](https://discourse.julialang.org/t/csc-kills-the-prospect-of-multithreading-shouldnt-julia-use-csr/102491) on the Julia Discourse.

To create matrices in the CSC format, the functions [csc_matrix](csc_matrix.md) and [csc_matrix_32](csc_matrix.md) can be used.
