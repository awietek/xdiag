---
title: Operator types
---

Generic operators in XDiag are represented as [OpSum](opsum.md) objects made up of a coupling, which can be a real/complex number or a string, and [Op](op.md) objects. Every [Op](op.md) is defined by a `type`. Here we list all the available types implemented in XDiag, their definition, their required number of sites, and the blocks for which they are available.

The precise definition of an operator can depend on the block it is applied to (for example, a `Hop` term sums over both spin species on a [tJ](../blocks/tJ.md) block, but acts on a single species on a [Fermion](../blocks/fermion.md) block). For the block-specific formulas, see the "Operators" section on the respective block page: [Spinhalf](../blocks/spinhalf.md#operators), [tJ](../blocks/tJ.md#operators), [Electron](../blocks/electron.md#operators), [Boson](../blocks/boson.md#operators), [Fermion](../blocks/fermion.md#operators).

Throughout, $c^\dagger_{i\sigma}, c_{i\sigma}$ denote fermionic creation/annihilation operators, $a^\dagger_i, a_i$ bosonic ones, $n_i$ a number operator, and $\mathbf{S}_i = (S^x_i, S^y_i, S^z_i)$ a spin operator with $S^\pm_i = S^x_i \pm i S^y_i$.

### Spin operators

| Type              | Description                                                              | No. of sites | Blocks                        |
|:------------------|:-------------------------------------------------------------------------|:------------:|:------------------------------|
| `Sz`              | local magnetic moment $S^z_i$                                            | 1            | Spinhalf, tJ, Electron, Boson |
| `Sx`              | local magnetic moment $S^x_i$                                            | 1            | Spinhalf, tJ, Electron, Boson |
| `Sy`              | local magnetic moment $S^y_i$                                            | 1            | Spinhalf, tJ, Electron, Boson |
| `S+`              | spin raising operator $S^+_i$                                            | 1            | Spinhalf, tJ, Electron, Boson |
| `S-`              | spin lowering operator $S^-_i$                                           | 1            | Spinhalf, tJ, Electron, Boson |
| `SdotS`           | Heisenberg interaction $\mathbf{S}_i \cdot \mathbf{S}_j$                 | 2            | Spinhalf, tJ, Electron, Boson |
| `SzSz`            | Ising interaction $S^z_i S^z_j$                                          | 2            | Spinhalf, tJ, Electron, Boson |
| `Exchange`        | spin exchange $\frac{1}{2}(S^+_i S^-_j +  S^-_i S^+_j)$                  | 2            | Spinhalf, tJ, Electron, Boson |
| `ExchangeAsym`    | antisymmetric exchange $\frac{1}{2}(S^+_i S^-_j - S^-_i S^+_j)$          | 2            | Spinhalf, tJ, Electron, Boson |
| `ScalarChirality` | scalar chirality $\mathbf{S}_i \cdot (\mathbf{S}_j \times \mathbf{S}_k)$ | 3            | Spinhalf, Boson               |
| `TotalSz`         | total magnetization $\sum_i S^z_i$                                       | 0            | Spinhalf, tJ, Electron        |

On a [Boson](../blocks/boson.md) block the spin operators refer to the spin $S = (d-1)/2$ associated with the local dimension $d$.

### Bosonic operators

| Type       | Description                                              | No. of sites | Blocks |
|:-----------|:---------------------------------------------------------|:------------:|:-------|
| `Adag`     | bosonic creation operator $a^\dagger_i$                  | 1            | Boson  |
| `A`        | bosonic annihilation operator $a_i$                      | 1            | Boson  |

### Fermionic operators

| Type      | Description                                              | No. of sites | Blocks                 |
|:----------|:---------------------------------------------------------|:------------:|:-----------------------|
| `Cdag`    | creation operator $c^\dagger_i$ (spinless)               | 1            | Fermion                |
| `C`       | annihilation operator $c_i$ (spinless)                   | 1            | Fermion                |
| `Cdagup`  | creation operator $c^\dagger_{i\uparrow}$                | 1            | tJ, Electron           |
| `Cup`     | annihilation operator $c_{i\uparrow}$                    | 1            | tJ, Electron           |
| `Cdagdn`  | creation operator $c^\dagger_{i\downarrow}$              | 1            | tJ, Electron           |
| `Cdn`     | annihilation operator $c_{i\downarrow}$                  | 1            | tJ, Electron           |

### Hopping operators

The hopping terms are the hermitian combination $-(c^\dagger_i c_j + c^\dagger_{j\sigma}c_{i\sigma})$ carrying an overall minus sign. Note that the coupling multiplies the whole term, so a complex coupling makes it non-hermitian (see [Complex couplings](opsum.md#complex-couplings)). The `*Asym` variants provide the antisymmetric combination $-(c^\dagger_i c_j - c^\dagger_{j\sigma}c_{i\sigma})$.

| Type        | Description                                                                                                              | No. of sites | Blocks                       |
|:------------|:-------------------------------------------------------------------------------------------------------------------------|:------------:|:-----------------------------|
| `Hop`       | hopping over all species, $-\sum_\sigma (c^\dagger_{i\sigma}c_{j\sigma} + c^\dagger_{j\sigma}c_{i\sigma})$               | 2            | tJ, Electron, Boson, Fermion |
| `Hopup`     | hopping of $\uparrow$ electrons, $-(c^\dagger_{i\uparrow}c_{j\uparrow} + c^\dagger_{j\uparrow}c_{i\uparrow})$            | 2            | tJ, Electron                 |
| `Hopdn`     | hopping of $\downarrow$ electrons, $-(c^\dagger_{i\downarrow}c_{j\downarrow} + c^\dagger_{j\downarrow}c_{i\downarrow})$  | 2            | tJ, Electron                 |
| `HopAsym`   | antisymmetric hopping over all species, $-\sum_\sigma (c^\dagger_{i\sigma}c_{j\sigma} - c^\dagger_{j\sigma}c_{i\sigma})$ | 2            | tJ, Electron, Boson, Fermion |
| `HopupAsym` | antisymmetric $\uparrow$ hopping, $-(c^\dagger_{i\uparrow}c_{j\uparrow} - c^\dagger_{j\uparrow}c_{i\uparrow})$           | 2            | tJ, Electron                 |
| `HopdnAsym` | antisymmetric $\downarrow$ hopping, $-(c^\dagger_{i\downarrow}c_{j\downarrow} - c^\dagger_{j\downarrow}c_{i\downarrow})$ | 2            | tJ, Electron                 |

### Density and interaction operators

| Type          | Description                                                        | No. of sites | Blocks                 |
|:--------------|:-------------------------------------------------------------------|:------------:|:-----------------------|
| `N`           | number operator $n_i$ (spinless fermion / boson)                   | 1            | Boson, Fermion         |
| `NN`          | density-density interaction $n_i n_j$ (spinless)                   | 2            | Fermion                |
| `Nup`         | number of $\uparrow$ electrons $n_{i\uparrow}$                     | 1            | tJ, Electron           |
| `Ndn`         | number of $\downarrow$ electrons $n_{i\downarrow}$                 | 1            | tJ, Electron           |
| `Ntot`        | total number $n_i = n_{i\uparrow} + n_{i\downarrow}$               | 1            | tJ, Electron           |
| `Nupdn`       | double occupancy $n_{i\uparrow} n_{i\downarrow}$                   | 1            | Electron               |
| `NtotNtot`    | density-density interaction $n_i n_j$                              | 2            | tJ, Electron           |
| `NupNup`      | $n_{i\uparrow} n_{j\uparrow}$                                      | 2            | tJ, Electron           |
| `NupNdn`      | $n_{i\uparrow} n_{j\downarrow}$                                    | 2            | tJ, Electron           |
| `NdnNup`      | $n_{i\downarrow} n_{j\uparrow}$                                    | 2            | tJ, Electron           |
| `NdnNdn`      | $n_{i\downarrow} n_{j\downarrow}$                                  | 2            | tJ, Electron           |
| `NupdnNupdn`  | double-occupancy correlation $n_{i\uparrow}n_{i\downarrow} n_{j\uparrow}n_{j\downarrow}$ | 2 | Electron   |
| `HubbardU`    | on-site interaction: $\sum_i n_{i\uparrow}n_{i\downarrow}$ (Electron), $\frac{1}{2}\sum_i n_i(n_i-1)$ (Boson) | 0 | Electron, Boson |

### $t$-$J$ operators

| Type       | Description                                                      | No. of sites | Blocks |
|:-----------|:-----------------------------------------------------------------|:------------:|:-------|
| `tJSzSz`   | $t$-$J$ Ising interaction $S^z_i S^z_j - \frac{n_i n_j}{4}$       | 2            | tJ     |
| `tJSdotS`  | $t$-$J$ Heisenberg interaction $\mathbf{S}_i \cdot \mathbf{S}_j - \frac{n_i n_j}{4}$ | 2 | tJ |

### Total quantum numbers

| Type        | Description                                     | No. of sites | Blocks                 |
|:------------|:-------------------------------------------------|:------------:|:-----------------------|
| `TotalN`    | total particle number $\sum_i n_i$               | 0            | tJ, Electron, Boson, Fermion |
| `TotalNup`  | total $\uparrow$ number $\sum_i n_{i\uparrow}$   | 0            | tJ, Electron           |
| `TotalNdn`  | total $\downarrow$ number $\sum_i n_{i\downarrow}$ | 0          | tJ, Electron           |

### Generic operators

| Type     | Description                                                       | No. of sites | Blocks                                |
|:---------|:------------------------------------------------------------------|:------------:|:--------------------------------------|
| `Id`     | the identity operator $\mathbb{1}$                                | 0            | Spinhalf, tJ, Electron, Boson, Fermion |
| `Matrix` | a generic interaction defined by an explicit matrix (see below)   | any          | Spinhalf, Boson                        |

!!! info "Distributed blocks"

	The distributed blocks ([SpinhalfDistributed](../blocks/spinhalf_distributed.md), [tJDistributed](../blocks/tJ_distributed.md), [ElectronDistributed](../blocks/electron_distributed.md)) support only the subset of operators that have a dedicated distributed kernel. The supported types are listed on each distributed block's page.

## Matrix type

The `Matrix` interaction type is a special type with which one can define generic interactions for the [Spinhalf](../blocks/spinhalf.md) and [Boson](../blocks/boson.md) blocks. In addition to the `type` and `sites` argument, a numerical matrix is provided when constructing the [Op](op.md) object. The matrix describes the operator acting on the $d^n$ dimensional space spanned by the $n$ sites of the operator (with local dimension $d$). For example, we can represent a $S^x$ spin operator as,

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:matrixtype1"
	```

More generically, we can use this mechanism to construct arbitrary spin interactions, e.g.

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:matrixtype2"
	```

Here we have been using the Kronecker product function `kron`.
