### Sparse Matrix Capabilities

Working with many-body quantum systems often involves matrices which only have a small number of non-zero elements, also known as *sparse matrices*.
Since storing and handling such objects in the conventional way (i.e. element by element) is very inefficient, XDiag includes basic implementations of three common [sparse-matrix types](documentation/algebra/sparse/sparse_matrix_types.md): the coordinate (COO), the compressed-sparse-row (CSR), and the compressed-sparse-column (CSC) formats.

Just as the [matrix](documentation/algebra/matrix.md) function can be used to obtain the full matrix representing a given operator `ops` (in the form of an [OpSum](documentation/operators/opsum.md)) on a given Hilbert space `block` (see section [Matrix representation](#matrix-representation)), there are functions [coo_matrix](documentation/algebra/sparse/coo_matrix.md), [csr_matrix](documentation/algebra/sparse/csr_matrix.md), and [csc_matrix](documentation/algebra/sparse/csc_matrix.md) to obtain the same matrix in the respective sparse format.

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:usage_guide_spm1"
	```
	
=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:usage_guide_spm1"
	```

Note that the C++ implementation distinguishes between real and complex matrices, e.g., there are the [csr_matrix](documentation/algebra/sparse/csr_matrix.md) and [csr_matrixC](documentation/algebra/sparse/csr_matrix.md) functions.

The objects returned by these functions are "raw" in the sense that they are not an instance of a sparse matrix implementation by another library, but contain all the information to call the respective sparse-matrix constructor of your sparse-matrix library of choice.
For instance, the output of [csc_matrix](documentation/algebra/sparse/csc_matrix.md) can be used to construct the `sp_mat` type implemented by the C++ [Armadillo library](https://arma.sourceforge.net) or the `SparseMatrixCSC` type implemented by [SparseArrays](https://docs.julialang.org/en/v1/stdlib/SparseArrays/) in julia.

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:usage_guide_spm2"
	```
	
=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:usage_guide_spm2"
	```

While the COO, CSR, and CSC formats are supported for extracting sparse matrices from XDiag, only the CSR format is used internally because it is the only one suitable for parallelized matrix-vector multiplications.
The following julia example finds the ground state of an open-boundary Heisenberg chain and compares the default (matrix-free, i.e., on-the-fly) implementation to first converting the Hamiltonian to the CSR format.

=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:usage_guide_spm3"
	```

It can be observed that (at medium block sizes) the runtime of both versions is comparable, while the matrix-free variant will always require substantially less memory.
Note that the computational effort and memory required to construct the CSR representation grows significantly with system size such that the matrix-free version will generally run faster on larger systems.
The CSR implementations of XDiag functions are thus only recommended when the CSR matrix can be pre-computed once and needs to be reused frequently.
