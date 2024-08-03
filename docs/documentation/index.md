---
title: Overview
---

# Documentation

## Algorithms


| Name                             | Description                                                         |                          Language |
|:---------------------------------|:--------------------------------------------------------------------|----------------------------------:|
| [eigval0](algorithms/eigval0.md) | Computes the lowest lying eigenvalue of an operator                 | :simple-cplusplus: :simple-julia: |
| [eig0](algorithms/eig0.md)       | Computes the lowest lying eigenvalue and eigenvector of an operator | :simple-cplusplus: :simple-julia: |

## Algebra
|                             |                                                                  |                                   |
|:----------------------------|:-----------------------------------------------------------------|----------------------------------:|
| [matrix](algebra/matrix.md) | Creates the full matrix representation of an operator on a block | :simple-cplusplus: :simple-julia: |

## Blocks
|                                |                                            |                                   |
|:-------------------------------|:-------------------------------------------|----------------------------------:|
| [Spinhalf](blocks/spinhalf.md) | Block of a spin $S=1/2$ type Hilbert space | :simple-cplusplus: :simple-julia: |
| [tJ](blocks/tJ.md)             | Block of a $t-J$ type Hilbert space        | :simple-cplusplus: :simple-julia: |
| [Electron](blocks/electron.md) | Block of a Electron type Hilbert space     | :simple-cplusplus: :simple-julia: |

## Operators
|                                   |                                                  |                                   |
|:----------------------------------|:-------------------------------------------------|----------------------------------:|
| [Op](operators/op.md)             | A local operator acting on several lattice sites | :simple-cplusplus: :simple-julia: |
| [OpSum](operators/opsum.md)       | Sum of local operators                           | :simple-cplusplus: :simple-julia: |
| [Coupling](operators/coupling.md) | Describes the coupling of a local operator       | :simple-cplusplus: :simple-julia: |


## States
|                          |                                                    |                                   |
|:-------------------------|:---------------------------------------------------|----------------------------------:|
| [State](states/state.md) | A generic state describing a quantum wave function | :simple-cplusplus: :simple-julia: |


## Symmetries
|                                                     |                                                     |                                   |
|:----------------------------------------------------|:----------------------------------------------------|----------------------------------:|
| [Permutation](symmetries/permutation.md)            | Permutations of indices or lattice sites            | :simple-cplusplus: :simple-julia: |
| [PermutationGroup](symmetries/permutation_group.md) | A group of permutations                             | :simple-cplusplus: :simple-julia: |
| [Representation](symmetries/representation.md)      | A (1D) irreducible representation of a finite group | :simple-cplusplus: :simple-julia: |

## Utilities

|                                       |                                               |                                   |
|:--------------------------------------|:----------------------------------------------|----------------------------------:|
| [Logging](utilities/logging.md)       | Controling what is written to standard output | :simple-cplusplus: :simple-julia: |
| [Timing](utilities/timing.md)         | Measurng wall time straightforwardly          |                :simple-cplusplus: |
| [XDIAG_SHOW](utilities/xdiag_show.md) | Macro for printing debugging information      |                :simple-cplusplus: |
