---
title: Releases
---

## v0.5.0

**Jul. 14, 2026**

Major release introducing a full operator algebra, two new Hilbert space types, a block eigensolver, and improved memory efficiency.

This release marks a significant step forward in both the expressiveness and the performance of XDiag. Operators are now first-class algebraic objects that can be multiplied and conjugated freely, two entirely new families of degrees of freedom — spinless fermions and bosons (including general spin-$S$) — are supported, and excited states can be computed reliably with a dedicated block eigensolver. At the same time, the memory footprint of every block type has been substantially reduced.

### New features

* **Arbitrary operators and a full operator algebra.**
  [OpSum](documentation/operators/opsum.md) objects can now not only be added, subtracted, and scaled, but also **multiplied** with one another, forming the (generally non-commutative) product of two operators. Together with the hermitian conjugation [hc](documentation/operators/hc.md), which acts as an involution, the OpSums now realize the mathematical structure of an **involutive algebra** ($*$-algebra). This makes it straightforward to build arbitrary composite operators — for example generic multi-point correlation functions such as $S^x_i S^y_j$ — directly from elementary building blocks. New convenience routines [expect](documentation/states/expect.md) and [correlation_matrix](documentation/states/correlation_matrix.md) evaluate local expectation values and two-point correlations across a whole lattice in a single call.

* **New Fermion and Boson blocks.**
  Two new Hilbert space types have been added: the [Fermion](documentation/blocks/fermion.md) block for spinless fermions, and the [Boson](documentation/blocks/boson.md) block for bosons with a configurable local dimension $d$. Since a $d$-level local degree of freedom equally describes a spin $S = (d-1)/2$, the Boson block (also available under the alias `Spin`) additionally provides general spin-$S$ physics, complete with the corresponding spin operators.

* **LOBPCG block eigensolver.**
  A new implementation of the **LOBPCG** algorithm (*Locally Optimal Block Preconditioned Conjugate Gradient*) has been added, exposed through the [eigvals](documentation/linalg/eigvals.md), [eigs](documentation/linalg/eigs.md), and [eigs_lobpcg](documentation/linalg/eigs_lobpcg.md) functions. As a block eigensolver, it reliably resolves several of the lowest eigenstates at once — including their degeneracies — which the Lanczos algorithm handles less robustly.

* **Improved memory efficiency throughout all blocks.**
  The internal representation of all block types has been reworked, substantially reducing the memory required to store and iterate over the basis states. This allows larger system sizes to be reached within the same memory budget.

### Breaking changes

* **Complex couplings on `Hop` and `Exchange` are no longer hermitian.**
  Previously, hopping (`Hop`, `Hopup`, `Hopdn`) and exchange (`Exchange`) terms were defined to remain hermitian even with a complex coupling. A complex coupling is now a plain prefactor that is complex-conjugated under [hc](documentation/operators/hc.md), which is required for the OpSums to form an algebra. To obtain the antisymmetric (non-hermitian) variants, the dedicated operator types `HopAsym`, `HopupAsym`, `HopdnAsym`, and `ExchangeAsym` have been introduced. See [Complex couplings](documentation/operators/opsum.md#complex-couplings) for details.

* **`ProductState` is now a vector of integers instead of strings.**
  The local configurations of a [ProductState](documentation/states/product_state.md) are now represented by integer local quantum numbers (e.g. `0` = $\downarrow$, `1` = $\uparrow$ for [Spinhalf](documentation/blocks/spinhalf.md)) rather than strings such as `"Up"`/`"Dn"`. Accordingly, `product_state` is constructed from an integer vector, and `to_string(pstate)` prints these integers; use `to_string(pstate, block)` to obtain the human-readable configuration. The integer encoding for each block type is documented on the respective block page.

## v0.4.1

**Dec. 4, 2025**
Small patch release

* Fixing bug reported in issue https://github.com/awietek/xdiag/issues/99

## v0.4.0

**Nov. 4, 2025**
Major update introducing sparse matrix capabilities

* Novel functionality for sparse matrices
* Sparse matrices can be created in COO, CSC, and CSR formats
* Internal iterative algorithms can be used with sparse matrices
* several minor bugfixes
* deprecated support for Julia versions 1.9 and lower


## v0.3.3

**Jun. 17, 2025**

Minor improvements and additional density operators

* Multithreaded parallelization works on aarch64 MacOS now (Julia version)
* introduced new density operators NupNup, NupNdn, NdnNup, NdnNdn for Electron block
* Exceptions thrown in constructors are now handed on to Julia
* OpSums can be created with integer coefficients
* several minor bugfixes

## v0.3.2

**Apr. 3, 2025**

ElectronDistributed and further enhancements

* Adds new block type: ElectronDistributed
* Added benchmarks
* linking to threaded Intel MKL enabled
* several enhancements and bugfixes

## v0.3.1

**Feb. 11, 2025**

New API and operator logic

* Updated to a new streamlined API
* Implemented operator logic including symmetries
* Completed wrapper for Julia
* w = apply(ops, v) works now

---

## v0.2.3

**Sep. 9, 2024**

Introduced 1-indexing everywhere in Julia version

* only changes to XDiag.jl, C++ untouched
* XDiag_jll.jl remains at v0.2.2

---

## v0.2.2

**Aug. 27, 2024**

Lanczos routines and multicolumn States

* wrapped eigs_lanczos, eigvals_lanczos
* implemented apply for States with multiple columns
* changed wrapping of std::vectors of Op and Permutation

---

## v0.2.1

**Aug. 16, 2024**

Small patch release providing small utility functions

* wrapped say_hello, print_version, and set_verbosity
* resorted to compiling wrapper with conventional OpenMP on aarch64 apple
* Fixed faulty behaviour of OpenMP on aarch64 apple

---

## v0.2.0

**Aug. 15, 2024**

Basic functionality for three Hilbert space types, Spinhalf, tJ, and Electron, has been implemented. Features are:

* Algebra with and without permutation symmetries
* Parallelization with OpenMP and MPI
* CMake has been properly set up
* Iterative algorithms present, Lanczos, Arnoldi, time evolution
* A minimal Julia wrapper has been written
* The Julia wrapper compiles on several target 64bit systems using BinaryBuilder
