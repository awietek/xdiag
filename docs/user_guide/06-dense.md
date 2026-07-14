---
title: Dense matrices
---

# Dense matrices

Given an operator in the form of an [OpSum](../documentation/operators/opsum.md) object and a Hilbert space (block), a dense matrix representation of the operator on the computational basis of the block can be computed using the [matrix](../documentation/kernels/matrix.md) function.

=== "Julia"
	```julia
	--8<-- "examples/user_guide/main.jl:usage_guide_mat1"
	```
=== "C++"
	```c++
	--8<-- "examples/user_guide/main.cpp:usage_guide_mat1"
	```

In C++, XDiag is using the [Armadillo library](https://arma.sourceforge.net/) with the `arma` namespace. The Armadillo library serves as the linear algebra backend of XDiag and can be used to perform further calculations. The [matrix](../documentation/kernels/matrix.md) function returns a real matrix (`arma::mat`). Whenever the operator has complex matrix elements, the complex-valued counterpart [matrixC](../documentation/kernels/matrix.md) must be used instead, which returns an `arma::cx_mat`. In Julia, this distinction is not necessary: `matrix` automatically returns a real or complex matrix depending on the operator.

To compute all eigenvalues and eigenvectors of a Hamiltonian, i.e. to perform a full exact diagonalization, standard linear algebra routines can be used.

=== "Julia"
	```julia
	--8<-- "examples/user_guide/main.jl:usage_guide_mat2"
	```
=== "C++"
	```c++
	--8<-- "examples/user_guide/main.cpp:usage_guide_mat2"
	```

In Julia, the `eigen` and `Symmetric` functions are part of the [Linear Algebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) standard library. In C++, `arma::eig_sym` is the corresponding [Armadillo](https://arma.sourceforge.net/) routine for the eigendecomposition of a symmetric (or hermitian) matrix.
