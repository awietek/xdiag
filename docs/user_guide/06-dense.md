### Matrix representation

Given an operator in the form of an [OpSum](documentation/operators/opsum.md) object and a Hilbert space (block) in the form of a certain block, a dense matrix representation of the operator on the computational basis of the block can be computed using the [matrix](documentation/algebra/matrix.md) function.

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:usage_guide_mat1"
	```
	
=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:usage_guide_mat1"
	```

In C++, XDiag is using the [Armadillo library](https://arma.sourceforge.net/) with the `arma` namespace. The Armadillo library serves as the linear algebra backend of XDiag and can be used to perform further calculations. To compute all eigenvalues and eigenvectors of a Hamiltonian, i.e. to perform a full exact diagonalization, standard linear algebra routines can be used.

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:usage_guide_mat2"
	```
	
=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:usage_guide_mat2"
	```

In Julia, the `eigen` and `Symmetric` functions are part of the [Linear Algebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) standard library.
