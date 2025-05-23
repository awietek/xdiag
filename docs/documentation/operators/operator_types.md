---
title: Operator types
---

## List of operator types

Generic operators in XDiag are represented as [OpSum](opsum.md) objects made up of a coupling, which can be a real/complex number or a string, and [Op](op.md) objects. Every [Op](op.md) is defined by a `type`. Here we list all the available types implemented in XDiag, their required number of sites, and the blocks for which they are available.

| Type              | Description                                                                                                                                                                | No. of sites | Blocks                                                                          |
|:------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------|---------------------------------------------------------------------------------|
| `Hop`             | A hopping term for $\uparrow$ and $\downarrow$ spins of the form $$ \textcolor{red}{-}\sum_{\sigma=\uparrow\downarrow} (tc^\dagger_{i\sigma}c_{j\sigma} + \textrm{h.c.})$$ | 2            | tJ, Electron, tJDistributed, ElectronDistributed                                |
| `Hopup`           | A hopping term for $\uparrow$ spins of the form $$ \textcolor{red}{-}(tc^\dagger_{i\uparrow}c_{j\uparrow} + \textrm{h.c.})$$                                               | 2            | tJ, Electron, tJDistributed, ElectronDistributed                                |
| `Hopdn`           | A hopping term for $\downarrow$ spins of the form $$ \textcolor{red}{-}(tc^\dagger_{i\downarrow}c_{j\downarrow} + \textrm{h.c.})$$                                         | 2            | tJ, Electron, tJDistributed, ElectronDistributed                                |
| `HubbardU`        | A uniform Hubbard interaction across the full lattice of the form $$ \sum_i n_{i\uparrow}n_{i\downarrow}$$                                                                 | 0            | Electron, ElectronDistributed                                                   |
| `Cdagup`          | A fermionic creation operator for an $\uparrow$ spin $c^\dagger_{i\uparrow}$                                                                                               | 1            | tJ, Electron, tJDistributed, ElectronDistributed                                |
| `Cdagdn`          | A fermionic creation operator for an $\downarrow$ spin $c^\dagger_{i\downarrow}$                                                                                           | 1            | tJ, Electron, tJDistributed, ElectronDistributed                                |
| `Cup`             | A fermionic annihilation operator for an $\uparrow$ spin $c_{i\uparrow}$                                                                                                   | 1            | tJ, Electron, tJDistributed, ElectronDistributed                                |
| `Cdn`             | A fermionic annihilation operator for an $\downarrow$ spin $c_{i\downarrow}$                                                                                               | 1            | tJ, Electron, tJDistributed, ElectronDistributed                                |
| `Nup`             | A number operator for an $\uparrow$ spin $n_{i\uparrow}$                                                                                                                   | 1            | tJ, Electron, tJDistributed, ElectronDistributed                                |
| `Ndn`             | A number operator for an $\downarrow$ spin $n_{i\downarrow}$                                                                                                               | 1            | tJ, Electron, tJDistributed, ElectronDistributed                                |
| `Ntot`            | A number operator $n_i = n_{i\uparrow} + n_{i\downarrow}$                                                                                                                  | 1            | tJ, Electron, tJDistributed, ElectronDistributed                                |
| `Nupdn`           | double occupancy $d_i = n_{i\uparrow} n_{i\downarrow}$                                                                                                                     | 1            | Electron, ElectronDistributed                                                   |
| `NupdnNupdn`      | double occupancy correlation $d_id_j$                                                                                                                                      | 2            | Electron, ElectronDistributed                                                   |
| `NtotNtot`        | A density-density interaction $n_i n_j$                                                                                                                                    | 2            | tJ, Electron, tJDistributed, ElectronDistributed                                |
| `SdotS`           | A Heisenberg interaction of the form $$ \mathbf{S}_i \cdot \mathbf{S}_j = S^x_iS^x_j + S^y_iS^y_j + S^z_iS^z_j$$                                                           | 2            | Spinhalf, tJ, Electron, SpinhalfDistributed, tJDistributed, ElectronDistributed |
| `SzSz`            | An Ising interaction of the form $ S^z_i S^z_j $                                                                                                                           | 2            | Spinhalf, tJ, Electron, SpinhalfDistributed, tJDistributed, ElectronDistributed |
| `Exchange`        | A spin exchange interaction of the form $$ \frac{1}{2}(JS^+_i S^-_j + J^*S^-_iS^+_j)$$                                                                                     | 2            | Spinhalf, tJ, Electron, SpinhalfDistributed, tJDistributed, ElectronDistributed |
| `Sz`              | A local magnetic moment in the $z$-direction $ S^z_i$                                                                                                                      | 1            | Spinhalf, tJ, Electron, SpinhalfDistributed, tJDistributed, ElectronDistributed |
| `S+`              | A local spin raising operator $S^+_i$                                                                                                                                      | 1            | Spinhalf, SpinhalfDistributed                                                   |
| `S-`              | A local spin lowering operator $S^-_i$                                                                                                                                     | 1            | Spinhalf, SpinhalfDistributed                                                   |
| `ScalarChirality` | A scalar chirality interaction of the form $$ \mathbf{S}_i \cdot ( \mathbf{S}_j \times  \mathbf{S}_k)$$                                                                    | 3            | Spinhalf                                                                        |
| `tJSzSz`          | An Ising interaction as encountered in the $t-J$ model of the form $$  S^z_i S^z_j - \frac{n_i n_j}{4}$$                                                                   | 2            | tJ, tJDistributed                                                               |
| `tJSdotS`         | An Heisenberg  interaction as encountered in the $t-J$ model of the form $$  \mathbf{S}_i \cdot \mathbf{S}_j - \frac{n_i n_j}{4}$$                                         | 2            | tJ, tJDistributed                                                               |
| `Matrix`          | A generic spin interaction no an arbitrary number of sites defined via a coupling matrix                                                                                   | arbitrary    | Spinhalf                                                                        |
| `NupNdn`          | A generic density-density interaction $n_{i\uparrow}n_{j\downarrow}$                                                                                                       | 2            | Electron                                                                        |
| `NdnNup`          | A generic density-density interaction $n_{i\downarrow}n_{j\uparrow}$                                                                                                       | 2            | Electron                                                                        |
| `NupNup`          | A generic density-density interaction for same-spin electrons $n_{i\uparrow}n_{i\uparrow}$                                                                                 | 2            | Electron                                                                        |
| `NdnNdn`          | A generic density-density interaction for same-spin electrons $n_{i\downarrow}n_{i\downarrow}$                                                                             | 2            | Electron                                                                        |

## Matrix type

The `Matrix` interaction type is a special type with whom one can define generic interactions for the [Spinhalf](../blocks/spinhalf.md) block. In addition to the `type` and `sites` argument, also a numerical matrix is provided when constructing the [Op](op.md) object. The matrix describes the operator acting on the $2^n$ dimensional space spanned by the $n$ sites of the operator. For example, we can represent a $S^x$ spin operator as,

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:matrixtype1"
	```
	
More generically, we can use this mechanism to construct arbitary spin interactions, e.g.	

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:matrixtype2"
	```

Here we have been using the Kronecker product function `kron`.
