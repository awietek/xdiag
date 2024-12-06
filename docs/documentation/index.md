---
title: Overview
---

# Documentation

## Building
| Name                                                   | Description                                                      |                          Language |
|:-------------------------------------------------------|:-----------------------------------------------------------------|----------------------------------:|
| [Compilation](compilation/advanced_compilation.md)     | advanced settings for compilation of the C++ library using CMake |                :simple-cplusplus: |
| [Documentation](compilation/building_documentation.md) | how to build and work on the documentation locally               |                 :simple-markdown: |
| [Julia Wrapper](compilation/julia_wrapper.md)          | how to build and develop the Julia wrapper locally               | :simple-cplusplus: :simple-julia: |
|                                                        |                                                                  |                                   |

## Algorithms


| Name                                             | Description                                                                                    |                          Language |
|:-------------------------------------------------|:-----------------------------------------------------------------------------------------------|----------------------------------:|
| [eigval0](algorithms/eigval0.md)                 | Computes the lowest lying eigenvalue of an operator                                            | :simple-cplusplus: :simple-julia: |
| [eig0](algorithms/eig0.md)                       | Computes the lowest lying eigenvalue and eigenvector of an operator                            | :simple-cplusplus: :simple-julia: |
| [eigvals_lanczos](algorithms/eigvals_lanczos.md) | Performs an iterative eigenvalue calculation using the Lanczos algorithm                       | :simple-cplusplus: :simple-julia: |
| [eigs_lanczos](algorithms/eigs_lanczos.md) | Performs an iterative eigenvalue calculation building eigenvectors using the Lanczos algorithm | :simple-cplusplus: :simple-julia: |

## Algebra
|                                       |                                                                     |                                   |
|:--------------------------------------|:--------------------------------------------------------------------|----------------------------------:|
| [matrix](algebra/matrix.md)           | Creates the full matrix representation of an operator on a block    | :simple-cplusplus: :simple-julia: |
| [apply](algebra/apply.md)             | Applies an operator to a state $\vert w \rangle = O \vert v\rangle$ | :simple-cplusplus: :simple-julia: |
| [norm](algebra/algebra.md#norm)       | Computes the 2-norm of a state                                      | :simple-cplusplus: :simple-julia: |
| [norm1](algebra/algebra.md#norm1)     | Computes the 1-norm of a state                                      | :simple-cplusplus: :simple-julia: |
| [norminf](algebra/algebra.md#norminf) | Computes the $\infty$-norm of a state                               | :simple-cplusplus: :simple-julia: |
| [dot](algebra/algebra.md#dot)         | Computes the dot product between two states                         | :simple-cplusplus: :simple-julia: |
| [inner](algebra/algebra.md#inner)     | Computes an expectation value $\langle v \vert O \vert v \rangle$   | :simple-cplusplus: :simple-julia: |

## Blocks
|                                |                                            |                                   |
|:-------------------------------|:-------------------------------------------|----------------------------------:|
| [Spinhalf](blocks/spinhalf.md) | Block of a spin $S=1/2$ type Hilbert space | :simple-cplusplus: :simple-julia: |
| [tJ](blocks/tJ.md)             | Block of a $t-J$ type Hilbert space        | :simple-cplusplus: :simple-julia: |
| [Electron](blocks/electron.md) | Block of a Electron type Hilbert space     | :simple-cplusplus: :simple-julia: |

## Operators
|                                       |                                                                           |                                   |
|:--------------------------------------|:--------------------------------------------------------------------------|----------------------------------:|
| [Op](operators/op.md)                 | A linear operator acting on the Hilbert space                             |                :simple-cplusplus: |
| [OpSum](operators/opsum.md)           | Sum of couplings times operators                                          |                :simple-cplusplus: |
| [isreal](operators/isreal.md)         | Returns whether an Op or OpSum is a real operator.                        |                :simple-cplusplus: |
| [hc](operators/hc.md)                 | Returns the hermitian conjugate of an Op or OpSum.                        |                :simple-cplusplus: |
| [symmetrize](operators/symmetrize.md) | Symmetrizes an operator with respect to a permutation symmetry group      | :simple-cplusplus: :simple-julia: |
| [Scalar](operators/scalar.md)         | A scalar number which can be either real or complex                       |                :simple-cplusplus: |
| [Coupling](operators/coupling.md)     | Describes a coupling associated with an operator, either string or scalar |                :simple-cplusplus: |


## States
|                                           |                                                                   |                                   |
|:------------------------------------------|:------------------------------------------------------------------|----------------------------------:|
| [State](states/state.md)                  | A generic state describing a quantum wave function                | :simple-cplusplus: :simple-julia: |
| [ProductState](states/product_state.md)   | A product state of local configurations                           | :simple-cplusplus: :simple-julia: |
| [RandomState](states/random_state.md)     | A random state with normal distributed coefficients               | :simple-cplusplus: :simple-julia: |
| [fill](states/fill.md)                    | Fill a state with a given model state                             | :simple-cplusplus: :simple-julia: |
| [product](states/create_state.md#product) | Creates a filled product state                                    | :simple-cplusplus: :simple-julia: |
| [rand](states/create_state.md#rand)       | Create a filled random state with normal distributed coefficients | :simple-cplusplus: :simple-julia: |
| [zeros](states/create_state.md#zeros)     | Create a filled state with all zero entries                       | :simple-cplusplus: :simple-julia: |
| [zero](states/create_state.md#zero)       | Set all coefficients of a given state to zero                     | :simple-cplusplus: :simple-julia: |


## Symmetries
|                                                     |                                                     |                                   |
|:----------------------------------------------------|:----------------------------------------------------|----------------------------------:|
| [Permutation](symmetries/permutation.md)            | Permutations of indices or lattice sites            | :simple-cplusplus: :simple-julia: |
| [PermutationGroup](symmetries/permutation_group.md) | A group of permutations                             | :simple-cplusplus: :simple-julia: |
| [Representation](symmetries/representation.md)      | A (1D) irreducible representation of a finite group | :simple-cplusplus: :simple-julia: |

## Utilities

|                                                   |                                                          |                                   |
|:--------------------------------------------------|:---------------------------------------------------------|----------------------------------:|
| [set_verbosity](utilities/utils.md#set_verbosity) | Sets how much information is printed during computations | :simple-cplusplus: :simple-julia: |
| [say_hello](utilities/utils.md#say_hello)         | Prints a nice welcome message with version number        | :simple-cplusplus: :simple-julia: |
| [print_version](utilities/utils.md#print_version) | Prints the plain version number                          | :simple-cplusplus: :simple-julia: |
| [Logging](utilities/logging.md)       | Controling what is written to standard output | :simple-cplusplus: :simple-julia: |
| [Timing](utilities/timing.md)         | Measurng wall time straightforwardly          |                :simple-cplusplus: |
| [XDIAG_SHOW](utilities/xdiag_show.md) | Macro for printing debugging information      |                :simple-cplusplus: |
