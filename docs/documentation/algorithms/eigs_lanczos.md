---
title: eigs_lanczos
---

Performs an iterative eigenvalue calculation building eigenvectors using the Lanczos algorithm. Returns the tridiagonal matrix, eigenvalues, number of iterations and the stopping criterion. The Lanczos interations are performed twice, where at the second run the eigenvectors are built.

**Sources** [eigs_lanczos.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/algorithms/lanczos/eigs_lanczos.hpp), [eigs_lanczos.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/algorithms/lanczos/eigs_lanczos.cpp)

---

## Definition

The Lanczos algorithm can be run in two distinct ways:

1. A random intial state $|\psi_0\rangle = |r\rangle$ with normal distributed entries is used.

	=== "C++"

		```c++
		eigs_lanczos_result_t
		eigs_lanczos(OpSum const &ops, Block const &block, int64_t neigvals = 1,
		             double precision = 1e-12, int64_t max_iterations = 1000,
                     double deflation_tol = 1e-7, int64_t random_seed = 42);
		```

	=== "Julia"
	
		``` julia
		function eigs_lanczos(
			ops::OpSum,
			block::Block;
			neigvals::Int64 = 1,
			precision::Float64 = 1e-12,
			max_iterations::Int64 = 1000,
			deflation_tol::Float64 = 1e-7,
			random_seed::Int64 = 42,
		)
		```

2. The initial state $|\psi_0\rangle$ is explicitly specified

	=== "C++"
		```c++
		eigs_lanczos_result_t 
		eigs_lanczos(OpSum const &ops, State const &state0, int64_t neigvals = 1,
                     double precision = 1e-12, int64_t max_iterations = 1000,
                     double deflation_tol = 1e-7);
		```

---

## Parameters

| Name           | Description                                                                       | Default |
|:---------------|:----------------------------------------------------------------------------------|---------|
| ops            | [OpSum](../operators/opsum.md) defining the bonds of the operator                 |         |
| block          | block on which the operator is defined                                            |         |
| neigvals       | number of eigenvalues to converge                                                 | 1       |
| precision      | accuracy of the computed ground state                                             | 1e-12   |
| max_iterations | maximum number of iterations                                                      | 1000    |
| deflation_tol  | tolerance for deflation, i.e. breakdown of Lanczos due to Krylow space exhaustion | 1e-7    |
| random_seed    | random seed for setting up the initial vector                                     | 42      |

---

## Returns

A struct with the following entries

| Entry        | Description                                                                                                           |
|:-------------|:----------------------------------------------------------------------------------------------------------------------|
| alphas       | diagonal elements of the tridiagonal matrix                                                                           |
| betas        | off-diagonal elements of the tridiagonal matrix                                                                       |
| eigenvalues  | the computed Ritz eigenvalues of the tridiagonal matrix                                                               |
| eigenvectors | [State](../states/state.md) of shape $D \times $`neigvals` holding all low-lying eigenvalues up to `neigvals` |
| niterations  | number of iterations performed                                                                                        |
| criterion    | string denoting the reason why the algorithm stopped                                                                  |

---

## Convergence criterion

The algorithm terminates if the $k$-th ($k$ is the argument `neigvals`) approximate eigenvalue changes only by a fraction smaller than $\epsilon$ ($k$ is the argument `precision`), i.e.

$$ (\tilde{e}_k^{(n)} - \tilde{e}_k^{(n-1)}) / \tilde{e}_k^{(n)} < \epsilon.$$

Here, $\tilde{e}_k^{(n)}$ denotes the Lanczos approximation to the $k$-th eigenvalue after $n$ iterations.

---

## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:eigs_lanczos"
	```

