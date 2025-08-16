---
title: Documentation
---

XDiag uses the C++ library [Armadillo](https://arma.sourceforge.net) as a linear algebra backend. Documentation for linear algebra operations can, therefore, be found in the [Armadillo Documentation](https://arma.sourceforge.net/docs.html).

## Algebra
| Name                                  | Description                                                               |                          Language |
|:--------------------------------------|:--------------------------------------------------------------------------|----------------------------------:|
| [matrix](algebra/matrix.md)           | Creates the full matrix representation of an operator on a block          | :simple-cplusplus: :simple-julia: |
| [apply](algebra/apply.md)             | Applies an operator to a state $\vert \phi \rangle = O \vert \psi\rangle$ | :simple-cplusplus: :simple-julia: |
| [dot](algebra/algebra.md#dot)         | Computes the dot product between two states                               | :simple-cplusplus: :simple-julia: |
| [inner](algebra/algebra.md#inner)     | Computes an expectation value $\langle \psi \vert O \vert \psi \rangle$   | :simple-cplusplus: :simple-julia: |
| [norm](algebra/algebra.md#norm)       | Computes the 2-norm of a state                                            | :simple-cplusplus: :simple-julia: |
| [norm1](algebra/algebra.md#norm1)     | Computes the 1-norm of a state                                            | :simple-cplusplus: :simple-julia: |
| [norminf](algebra/algebra.md#norminf) | Computes the $\infty$-norm of a state                                     | :simple-cplusplus: :simple-julia: |

#### Sparse matrices
| Name                                                         | Description                                                                      |                          Language |
|:-------------------------------------------------------------|:---------------------------------------------------------------------------------|----------------------------------:|
| [Sparse matrix types](algebra/sparse/sparse_matrix_types.md) | Explanation of different sparse matrix types and their use                       | :simple-cplusplus: :simple-julia: |
| [coo_matrix](algebra/sparse/coo_matrix.md)                   | Creates sparse matrix an operator in the coordinate (COO) format                 | :simple-cplusplus: :simple-julia: |
| [csr_matrix](algebra/sparse/csr_matrix.md)                   | Creates sparse matrix an operator in the compressed-sparse-row (CSR) format      | :simple-cplusplus: :simple-julia: |
| [csc_matrix](algebra/sparse/csc_matrix.md)                   | Creates sparse matrix an operator in the compressed-sparse-column (CSC) format   | :simple-cplusplus: :simple-julia: |
| [apply](algebra/sparse/apply.md)                             | Sparse matrix-vector (and sparse matrix-matrix) multiplication with CSR matrices | :simple-cplusplus: :simple-julia: |

---

## Algorithms

#### Diagonalization

| Name                                             | Description                                                                                    |                          Language |
|:-------------------------------------------------|:-----------------------------------------------------------------------------------------------|----------------------------------:|
| [eigval0](algorithms/eigval0.md)                 | Computes the lowest lying eigenvalue of an operator                                            | :simple-cplusplus: :simple-julia: |
| [eig0](algorithms/eig0.md)                       | Computes the lowest lying eigenvalue and eigenvector of an operator                            | :simple-cplusplus: :simple-julia: |
| [eigvals_lanczos](algorithms/eigvals_lanczos.md) | Performs an iterative eigenvalue calculation using the Lanczos algorithm                       | :simple-cplusplus: :simple-julia: |
| [eigs_lanczos](algorithms/eigs_lanczos.md)       | Performs an iterative eigenvalue calculation building eigenvectors using the Lanczos algorithm | :simple-cplusplus: :simple-julia: |

#### Time evolution

| Name                                                         | Description                                                                                                                                     |                          Language |
|:-------------------------------------------------------------|:------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------:|
| [time_evolve](algorithms/time_evolve.md)                     | Performs a real-time evolution $e^{ -iHt} \vert \psi \rangle$ of a State with a given Hermitian operator $H$                                | :simple-cplusplus: :simple-julia: |
| [imaginary_time_evolve](algorithms/imaginary_time_evolve.md) | Performs a imaginary-time evolution $e^{ -\tau H}\vert\psi\rangle$ of a State with a given Hermitian operator $H$                               | :simple-cplusplus: :simple-julia: |
| [evolve_lanczos](algorithms/evolve_lanczos.md)               | Computes the exponential $e^{z H}\vert\psi\rangle $ of a Hermitian operator times a State for a real or complex $z$ using the Lanczos algorithm | :simple-cplusplus: :simple-julia: |
| [time_evolve_expokit](algorithms/time_evolve_expokit.md)     | Performs a real-time evolution $e^{ -iHt} \vert \psi \rangle$ using a highly accurate Lanczos algorithm                                     | :simple-cplusplus: :simple-julia: |


---

## Blocks

**Shared memory**

| Name                                                  | Description                                                        |                          Language |
|:------------------------------------------------------|:-------------------------------------------------------------------|----------------------------------:|
| [Spinhalf](blocks/spinhalf.md)                        | Block of a spin $S=1/2$ type Hilbert space                         | :simple-cplusplus: :simple-julia: |
| [tJ](blocks/tJ.md)                                    | Block of a $t-J$ type Hilbert space                                | :simple-cplusplus: :simple-julia: |
| [Electron](blocks/electron.md)                        | Block of a Electron type Hilbert space                             | :simple-cplusplus: :simple-julia: |

**Distributed memory**

| Name                                                  | Description                                                        |           Language |
|:------------------------------------------------------|:-------------------------------------------------------------------|-------------------:|
| [SpinhalfDistributed](blocks/spinhalf_distributed.md) | Block of a spin $S=1/2$ type Hilbert space (distributed computing) | :simple-cplusplus: |
| [tJDistributed](blocks/tJ_distributed.md)             | Block of a $t-J$ type Hilbert space  (distributed computing)       | :simple-cplusplus: |
| [ElectronDistributed](blocks/electron_distributed.md) | Block of a Electron type Hilbert space  (distributed computing)    | :simple-cplusplus: |

---

## Building
| Name                                                   | Description                                                      |                          Language |
|:-------------------------------------------------------|:-----------------------------------------------------------------|----------------------------------:|
| [Compilation](compilation/compilation.md)     | Advanced settings for compilation of the C++ library using CMake |                :simple-cplusplus: |
| [Documentation](compilation/building_documentation.md) | How to build and work on the documentation locally               |                 :simple-markdown: |
| [Julia Wrapper](compilation/julia_wrapper.md)          | How to build and develop the Julia wrapper locally               | :simple-cplusplus: :simple-julia: |

---

## Input / Output

**TOML**

| Name                                                   | Description                                                                  |                          Language |
|:-------------------------------------------------------|:-----------------------------------------------------------------------------|----------------------------------:|
| [FileToml](io/file_toml.md)                            | A file handler for TOML files                                                | :simple-cplusplus: :simple-julia: |
| [read_opsum](io/read_opsum.md)                         | reads an [OpSum](operators/opsum.md) from a TOML file                        | :simple-cplusplus: :simple-julia: |
| [read_permutation_group](io/read_permutation_group.md) | reads a [PermutationGroup](symmetries/permutation_group.md) from a TOML file | :simple-cplusplus: :simple-julia: |
| [read_representation](io/read_representation.md) | reads a [Representation](symmetries/representation.md) from a TOML file | :simple-cplusplus: :simple-julia: |

**HDF5**

| Name                    | Description                                                               |           Language |
|:------------------------|:--------------------------------------------------------------------------|-------------------:|
| [FileH5](io/file_h5.md) | A file handler for [hdf5](https://www.hdfgroup.org/solutions/hdf5/) files | :simple-cplusplus: |

---


## Operators
|                                               |                                                      |   |
|:----------------------------------------------|:-----------------------------------------------------|--:|
| [Operator types](operators/operator_types.md) | A summary of all the operator types defined in XDiag |   |

| Name                                  | Description                                                                                                                          |                           Language |
|:--------------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------:|
| [Op](operators/op.md)                 | A linear operator acting on the Hilbert space                                                                                        | :simple-cplusplus:  :simple-julia: |
| [OpSum](operators/opsum.md)           | Sum of couplings times operators                                                                                                     | :simple-cplusplus:  :simple-julia: |
| [hc](operators/hc.md)                 | Returns the hermitian conjugate of an Op or OpSum.                                                                                   | :simple-cplusplus:  :simple-julia: |
| [symmetrize](operators/symmetrize.md) | Symmetrizes an operator with a [PermutationGroup](symmetries/permutation_group.md) or [Representation](symmetries/representation.md) |  :simple-cplusplus: :simple-julia: |

---

## States
| Name                                                  | Description                                                       |                          Language |
|:------------------------------------------------------|:------------------------------------------------------------------|----------------------------------:|
| [State](states/state.md)                              | A generic state describing a quantum wave function                | :simple-cplusplus: :simple-julia: |
| [ProductState](states/product_state.md)               | A product state of local configurations                           | :simple-cplusplus: :simple-julia: |
| [RandomState](states/random_state.md)                 | A random state with normal distributed coefficients               | :simple-cplusplus: :simple-julia: |
| [fill](states/fill.md)                                | Fill a state with a given model state                             | :simple-cplusplus: :simple-julia: |
| [product_state](states/create_state.md#product_state) | Creates a filled product state                                    | :simple-cplusplus: :simple-julia: |
| [random_state](states/create_state.md#random_state)   | Create a filled random state with normal distributed coefficients | :simple-cplusplus: :simple-julia: |
| [zero_state](states/create_state.md#zero_state)       | Create a filled state with all zero entries                       | :simple-cplusplus: :simple-julia: |
| [zero](states/create_state.md#zero)                   | Set all coefficients of a given state to zero                     | :simple-cplusplus: :simple-julia: |

---

## Symmetries
| Name                                                | Description                                         |                          Language |
|:----------------------------------------------------|:----------------------------------------------------|----------------------------------:|
| [Permutation](symmetries/permutation.md)            | Permutations of indices or lattice sites            | :simple-cplusplus: :simple-julia: |
| [PermutationGroup](symmetries/permutation_group.md) | A group of permutations                             | :simple-cplusplus: :simple-julia: |
| [Representation](symmetries/representation.md)      | A (1D) irreducible representation of a finite group | :simple-cplusplus: :simple-julia: |

---

## Utilities

| Name                                              | Description                                              |                          Language |
|:--------------------------------------------------|:---------------------------------------------------------|----------------------------------:|
| [set_verbosity](utilities/utils.md#set_verbosity) | Sets how much information is printed during computations | :simple-cplusplus: :simple-julia: |
| [say_hello](utilities/utils.md#say_hello)         | Prints a nice welcome message with version number        | :simple-cplusplus: :simple-julia: |
| [print_version](utilities/utils.md#print_version) | Prints the plain version number                          | :simple-cplusplus: :simple-julia: |
| [Logging](utilities/logging.md)                   | Controling what is written to standard output            | :simple-cplusplus: :simple-julia: |
| [Timing](utilities/timing.md)                     | Measurng wall time straightforwardly                     |                :simple-cplusplus: |
| [XDIAG_SHOW](utilities/xdiag_show.md)             | Macro for printing debugging information                 |                :simple-cplusplus: |
