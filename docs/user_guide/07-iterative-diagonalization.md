### Iterative algorithms

XDiag features several iterative algorithms to perform diagonalization and time evolution, which do not require dense matrix storage of the involved operators. Instead, applications of operators are implemented *on-the-fly* (i.e. matrix-free) and to minimize memory requirements.

#### Diagonalization

A fundamental property of a quantum system is its ground state energy, which in XDiag can be easily computed using the [eigval0](documentation/algorithms/eigval0.md) function.

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:usage_guide_iter1"
	```
	
=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:usage_guide_iter1"
	```

Similarly, the ground state can be computed using the [eig0](documentation/algorithms/eig0.md) function.

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:usage_guide_iter2"
	```
	
=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:usage_guide_iter2"
	```

Here, `e0` is a double precision real number and `psi0` is a [State](documentation/states/state.md) object. We employ the [Lanczos algorithm](https://nvlpubs.nist.gov/nistpubs/jres/045/jresv45n4p255_A1b.pdf) to perform iterative diagonalizations. While the functions [eigval0](documentation/algorithms/eigval0.md) and [eig0](documentation/algorithms/eig0.md) are convenient, XDiag offers also more refined routines called [eigvals_lanczos](documentation/algorithms/eigvals_lanczos.md) and [eigs_lanczos](documentation/algorithms/eigs_lanczos.md) which can target higher excited states and offer more control over the convergence properties. Moreover, they return the convergence criterion as well as the tridiagonal matrix in the Lanczos algorithm, which contains more information than only extremal eigenvalues. The functions [eig0](documentation/algorithms/eig0.md) and [eigs_lanczos](documentation/algorithms/eigs_lanczos.md) for computing eigenvalues perform the Lanczos iteration twice, first to compute the tridiagonal matrix and in a second run to build the corresponding eigenvectors in order to minimize memory requirements. For a precise description of these methods, we refer to their respective documentations.
