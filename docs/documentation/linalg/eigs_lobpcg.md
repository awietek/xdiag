---
title: eigs_lobpcg
---

Computes the `neigs` algebraically smallest eigenvalues and eigenvectors of a hermitian operator with the **LOBPCG** algorithm (*Locally Optimal Block Preconditioned Conjugate Gradient*). LOBPCG is a *block* eigensolver: it iterates a whole set of trial vectors simultaneously and is therefore well suited to compute several of the lowest eigenpairs at once and to reliably resolve degeneracies. This function exposes the full result of the algorithm, including residual norms and convergence histories; the convenience wrappers [eigvals](eigvals.md) and [eigs](eigs.md) are built on top of it.

The algorithm can be run either *on-the-fly* (matrix-free) or using a *sparse matrix* in the compressed-sparse-row format (see [CSRMatrix](../kernels/sparse/sparse_matrix_types.md)).

**Sources:** [eigs_lobpcg.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/linalg/lobpcg/eigs_lobpcg.hpp) · [eigs_lobpcg.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/linalg/lobpcg/eigs_lobpcg.cpp)

---

## Definition

#### On-the-fly

=== "C++"
	```c++
	EigsLobpcgResult eigs_lobpcg(OpSum const &ops, Block const &block,
	                             int64_t neigs = 1, int64_t guard = 2,
	                             double tol = 1e-10, int64_t max_iterations = 1000,
	                             int64_t random_seed = 42);
	```
=== "Julia"
	```julia
	eigs_lobpcg(ops::OpSum, block::Block; neigs::Int64 = 1, guard::Int64 = 2,
	            tol::Float64 = 1e-10, max_iterations::Int64 = 1000,
	            random_seed::Int64 = 42)
	```

#### Sparse matrix

=== "C++"
	```c++
	template <typename idx_t, typename coeff_t>
	EigsLobpcgResult eigs_lobpcg(CSRMatrix<idx_t, coeff_t> const &A,
	                             Block const &block, int64_t neigs = 1,
	                             int64_t guard = 2, double tol = 1e-10,
	                             int64_t max_iterations = 1000,
	                             int64_t random_seed = 42);
	```
=== "Julia"
	```julia
	eigs_lobpcg(A::CSRMatrix, block::Block; neigs::Int64 = 1, guard::Int64 = 2,
	            tol::Float64 = 1e-10, max_iterations::Int64 = 1000,
	            random_seed::Int64 = 42)
	```

---

## Parameters

| Name           | Description                                                                                                                | Default |
|:---------------|:---------------------------------------------------------------------------------------------------------------------------|---------|
| ops / A        | [OpSum](../operators/opsum.md) or [CSRMatrix](../kernels/sparse/sparse_matrix_types.md) defining the bonds of the operator |         |
| block          | block on which the operator is defined                                                                                     |         |
| neigs          | number of (lowest) eigenpairs to compute                                                                                   | 1       |
| guard          | number of additional guard vectors iterated on top of `neigs`, so a degenerate multiplet sitting exactly at the `neigs`-th eigenvalue is captured with the correct multiplicity | 2 |
| tol            | convergence tolerance on the residual norms                                                                                | 1e-10   |
| max_iterations | maximum number of iterations                                                                                               | 1000    |
| random_seed    | random seed for the initial block of vectors                                                                              | 42      |

---

## Returns

A struct with the following entries

| Entry                  | Description                                                                                                          |
|:-----------------------|:---------------------------------------------------------------------------------------------------------------------|
| eigenvalues            | the `neigs` lowest eigenvalues in ascending order                                                                    |
| eigenvectors           | [State](../states/state.md) of shape $D \times$ `neigs` holding the corresponding eigenvectors                       |
| residual_norms         | the final residual norm $\Vert H|\psi\rangle - \varepsilon |\psi\rangle \Vert$ of every computed eigenvector         |
| niterations            | number of iterations performed                                                                                       |
| criterion              | string denoting the reason why the algorithm stopped                                                                 |
| eigenvalue_history     | the eigenvalue estimates recorded at every iteration (useful for monitoring convergence)                            |
| residual_norms_history | the residual norms recorded at every iteration                                                                       |

---

## Usage Example

=== "Julia"
	```julia
	block = Spinhalf(12)
	ops = OpSum()
	for i in 1:12
	    ops += "J" * Op("SdotS", [i, mod1(i + 1, 12)])
	end
	ops["J"] = 1.0
	res = eigs_lobpcg(ops, block; neigs = 3)
	@show res.eigenvalues
	@show res.residual_norms
	```
=== "C++"
	```c++
	auto block = Spinhalf(12);
	auto ops = OpSum();
	for (int i = 0; i < 12; ++i) {
	  ops += "J" * Op("SdotS", {i, (i + 1) % 12});
	}
	ops["J"] = 1.0;
	auto res = eigs_lobpcg(ops, block, 3);
	XDIAG_SHOW(res.eigenvalues);
	XDIAG_SHOW(res.residual_norms);
	```
