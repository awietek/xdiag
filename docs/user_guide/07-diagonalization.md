---
title: Diagonalization
---

# Diagonalization

XDiag features several iterative algorithms to perform diagonalization and time evolution, which do not require dense matrix storage of the involved operators. Instead, applications of operators are implemented *on-the-fly* (i.e. matrix-free) to minimize memory requirements. All routines described here therefore take the operator directly as an [OpSum](../documentation/operators/opsum.md) together with the block it acts on, rather than a precomputed matrix.

## Ground state

A fundamental property of a quantum system is its ground state energy, which in XDiag can be easily computed using the [eigval0](../documentation/linalg/eigval0.md) function.

=== "Julia"
	```julia
	--8<-- "examples/user_guide/main.jl:usage_guide_iter1"
	```
=== "C++"
	```c++
	--8<-- "examples/user_guide/main.cpp:usage_guide_iter1"
	```

Similarly, the ground state itself can be computed using the [eig0](../documentation/linalg/eig0.md) function.

=== "Julia"
	```julia
	--8<-- "examples/user_guide/main.jl:usage_guide_iter2"
	```
=== "C++"
	```c++
	--8<-- "examples/user_guide/main.cpp:usage_guide_iter2"
	```

Here, `e0` is a double precision real number and `psi0` is a [State](../documentation/states/state.md) object. These two functions are convenient shortcuts for the underlying **Lanczos** implementation described next: they set up a Lanczos iteration, converge to the extremal (lowest) eigenvalue, and hide the additional bookkeeping.

## The Lanczos interface

The [Lanczos algorithm](https://nvlpubs.nist.gov/nistpubs/jres/045/jresv45n4p255_A1b.pdf) builds up a Krylov space and projects the operator onto a small tridiagonal matrix, whose extremal eigenvalues rapidly approximate those of the full operator. XDiag exposes this machinery directly through the [eigvals_lanczos](../documentation/linalg/eigvals_lanczos.md) and [eigs_lanczos](../documentation/linalg/eigs_lanczos.md) functions, the latter also returning the eigenvectors.

=== "Julia"
	```julia
	--8<-- "examples/user_guide/main.jl:usage_guide_iter_lanczos"
	```
=== "C++"
	```c++
	--8<-- "examples/user_guide/main.cpp:usage_guide_iter_lanczos"
	```

These routines return a result object containing more than just the eigenvalues: the coefficients `alphas` and `betas` of the tridiagonal Lanczos matrix, the number of iterations `niterations`, and the convergence `criterion` that was met. They also offer more control over the convergence properties, for example through the `precision`, `max_iterations`, and `deflation_tol` parameters, and they can be started from a user-provided initial [State](../documentation/states/state.md) instead of a random one. The [eig0](../documentation/linalg/eig0.md) and [eigs_lanczos](../documentation/linalg/eigs_lanczos.md) functions perform the Lanczos iteration twice: first to compute the tridiagonal matrix, and in a second run to build the eigenvectors, in order to minimize memory requirements.

A limitation of the Lanczos algorithm is that, converging essentially one vector at a time, it is not well suited to resolve degeneracies: eigenvalues with a multiplicity greater than one typically show up only once, and nearly degenerate excited states can be difficult to separate reliably.

## Excited states

To reliably compute several of the lowest eigenstates, including their degeneracies, XDiag provides the [eigvals](../documentation/linalg/eigvals.md) and [eigs](../documentation/linalg/eigs.md) functions. [eigvals](../documentation/linalg/eigvals.md) computes the `neigs` algebraically smallest eigenvalues, and [eigs](../documentation/linalg/eigs.md) additionally returns the corresponding eigenvectors.

=== "Julia"
	```julia
	--8<-- "examples/user_guide/main.jl:usage_guide_iter_eigvals"
	```
=== "C++"
	```c++
	--8<-- "examples/user_guide/main.cpp:usage_guide_iter_eigvals"
	```

---

=== "Julia"
	```julia
	--8<-- "examples/user_guide/main.jl:usage_guide_iter_eigs"
	```
=== "C++"
	```c++
	--8<-- "examples/user_guide/main.cpp:usage_guide_iter_eigs"
	```

Here `eigenvalues` is an `arma::vec` (Julia: `Vector{Float64}`) holding the `neigs` lowest eigenvalues in ascending order, while `eigenvectors` is a single [State](../documentation/states/state.md) object holding the `neigs` eigenvectors as its columns.

These functions are based on the **LOBPCG** algorithm (*Locally Optimal Block Preconditioned Conjugate Gradient*). In contrast to Lanczos, LOBPCG is a **block** eigensolver that iterates a whole set of trial vectors simultaneously. This requires more memory, since several vectors have to be kept in storage at once, but in return it reliably resolves excited states and, in particular, their degeneracies.

## The LOBPCG interface

For finer control, the [eigs_lobpcg](../documentation/linalg/eigs_lobpcg.md) function exposes the algorithm directly.

=== "Julia"
	```julia
	--8<-- "examples/user_guide/main.jl:usage_guide_iter_lobpcg"
	```
=== "C++"
	```c++
	--8<-- "examples/user_guide/main.cpp:usage_guide_iter_lobpcg"
	```

It returns a result object collecting, besides the `eigenvalues` and `eigenvectors`, additional information about the run: the final `residual_norms` of every eigenvector, the number of iterations `niterations`, the convergence `criterion` that was met, and the full `eigenvalue_history` and `residual_norms_history` across iterations, which are useful for monitoring convergence. Beyond `neigs`, the tolerance `tol` and the maximal number of iterations `max_iterations` can be set. The `guard` parameter enlarges the iterated block to `neigs + guard` vectors, so that a degenerate multiplet sitting exactly at the `neigs`-th eigenvalue is still captured with the correct multiplicity.
