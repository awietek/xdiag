---
title: eigs_lanczos
---

Performs an iterative eigenvalue calculation building eigenvectors using the Lanczos algorithm 

**Source** [eigs_lanczos.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/algorithms/lanczos/eigs_lanczos.hpp)

=== "Julia"

	``` julia
	function eigs_lanczos(
		ops::OpSum,
		block::Block;
		neigvals::Int64 = 1,
		precision::Float64 = 1e-12,
		max_iterations::Int64 = 1000,
		force_complex::Bool = false,
		deflation_tol::Float64 = 1e-7,
		random_seed::Int64 = 42,
	)
	```

=== "C++"

    ```c++
    eigs_lanczos_result_t
	eigs_lanczos(OpSum const &ops, Block const &block, int64_t neigvals = 1,
	            double precision = 1e-12, int64_t max_iterations = 1000,
                bool force_complex = false, double deflation_tol = 1e-7,
                int64_t random_seed = 42);
	```

## Parameters

| Name           | Description                                                                       | Default |
|:---------------|:----------------------------------------------------------------------------------|---------|
| ops            | [OpSum](../operators/opsum.md) defining the bonds of the operator                 |         |
| block          | block on which the operator is defined                                            |         |
| neigvals       | number of eigenvalues to converge                                                 | 1       |
| precision      | accuracy of the computed ground state                                             | 1e-12   |
| max_iterations | maximum number of iterations                                                      | 1000    |
| force_complex  | whether or not computation should be forced to have complex arithmetic            | false   |
| deflation_tol  | tolerance for deflation, i.e. breakdown of Lanczos due to Krylow space exhaustion | 1e-7    |
| random_seed    | random seed for setting up the initial vector                                     | 42      |

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
	
