---
title: Electron
---

A block in an Electron (fermions with $\uparrow, \downarrow$ spin) Hilbert space. 

**Sources:** [electron.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/blocks/electron.hpp) · [electron.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/blocks/electron.cpp)

## Constructors

=== "C++"	
	```c++
    Electron(int64_t nsites, std::string backend = "auto");
    Electron(int64_t nsites, int64_t nup, int64_t ndn, std::string backend = "auto");
    Electron(int64_t nsites, Representation irrep, std::string backend = "auto");
    Electron(int64_t nsites, int64_t nup, int64_t ndn, Representation irrep, std::string backend = "auto");
	```
	
=== "Julia"
	```julia
	Electron(nsites::Int64, backend::String="auto")
	Electron(nsites::Int64, nup::Int64, ndn::Int64, backend::String="auto")
	Electron(nsites::Int64, irrep::Representation, backend::String="auto")
	Electron(nsites::Int64, nup::Int64, ndn::Int64, irrep::Representation, backend::String="auto")
	```

| Name    | Description                                                                          |        |
|:--------|:-------------------------------------------------------------------------------------|--------|
| nsites  | number of sites (integer)                                                            |        |
| nup     | number of "up" electrons (integer)                                                   |        |
| ndn     | number of "dn" electrons (integer)                                                   |        |
| irrep   | Irreducible [Representation](../symmetries/representation.md)  of the symmetry group |        |
| backend | backend used for coding the basis states                                             | `auto` |

The parameter `backend` chooses how the block is coded internally. By using the default parameter `auto` the backend is chosen automatically. Alternatives are `32bit`, `64bit`

## Local configurations

Each site of an Electron block carries a local dimension $d=4$, allowing double occupancy. In a [ProductState](../states/product_state.md), the local configuration of every site is given by an integer with the following meaning:

| Integer | Configuration          | Symbol |
|:--------|:-----------------------|:-------|
| `0`     | empty                  | ○      |
| `1`     | up-spin electron       | ↑      |
| `2`     | down-spin electron     | ↓      |
| `3`     | doubly occupied        | ⇅      |

The integer encodes the occupation bit-wise: bit $0$ is the $\uparrow$ occupation, bit $1$ the $\downarrow$ occupation. The Jordan-Wigner sign convention used for the fermionic operators is described in the [Hilbert spaces](../../user_guide/03-hilbert-spaces.md#normal-ordering-of-fermionic-blocks) section of the user guide.

## Operators

The following operator types can be used on an Electron block. Here $c^\dagger_{i\sigma}$, $c_{i\sigma}$ are the electron creation and annihilation operators, $n_{i\sigma} = c^\dagger_{i\sigma}c_{i\sigma}$, $n_i = n_{i\uparrow} + n_{i\downarrow}$, and $\mathbf{S}_i$ is the spin operator with $S^z_i = \tfrac12(n_{i\uparrow} - n_{i\downarrow})$, $S^+_i = c^\dagger_{i\uparrow}c_{i\downarrow}$.

| Type          | Description                              | Formula                                                                       | No. of sites |
|:--------------|:-----------------------------------------|:------------------------------------------------------------------------------|:------------:|
| `Cdagup`      | creation operator ($\uparrow$)           | $c^\dagger_{i\uparrow}$                                                       | 1            |
| `Cup`         | annihilation operator ($\uparrow$)       | $c_{i\uparrow}$                                                               | 1            |
| `Cdagdn`      | creation operator ($\downarrow$)         | $c^\dagger_{i\downarrow}$                                                     | 1            |
| `Cdn`         | annihilation operator ($\downarrow$)     | $c_{i\downarrow}$                                                             | 1            |
| `Hop`         | hopping ($\uparrow$ and $\downarrow$)    | $-\sum_\sigma (c^\dagger_{i\sigma}c_{j\sigma} + \mathrm{h.c.})$               | 2            |
| `Hopup`       | hopping ($\uparrow$)                     | $-(c^\dagger_{i\uparrow}c_{j\uparrow} + \mathrm{h.c.})$                       | 2            |
| `Hopdn`       | hopping ($\downarrow$)                   | $-(c^\dagger_{i\downarrow}c_{j\downarrow} + \mathrm{h.c.})$                   | 2            |
| `HopAsym`, `HopupAsym`, `HopdnAsym` | antisymmetric hopping variants | $-(c^\dagger_{i\sigma}c_{j\sigma} - \mathrm{h.c.})$              | 2            |
| `HubbardU`    | on-site Hubbard interaction              | $\sum_i n_{i\uparrow} n_{i\downarrow}$                                        | 0            |
| `Nup`         | number operator ($\uparrow$)             | $n_{i\uparrow}$                                                              | 1            |
| `Ndn`         | number operator ($\downarrow$)           | $n_{i\downarrow}$                                                            | 1            |
| `Ntot`        | number operator                          | $n_i = n_{i\uparrow} + n_{i\downarrow}$                                       | 1            |
| `Nupdn`       | double occupancy                         | $n_{i\uparrow} n_{i\downarrow}$                                              | 1            |
| `NtotNtot`    | density-density interaction              | $n_i n_j$                                                                     | 2            |
| `NupdnNupdn`  | double-occupancy correlation             | $n_{i\uparrow}n_{i\downarrow}\, n_{j\uparrow}n_{j\downarrow}$                 | 2            |
| `NupNup`, `NupNdn`, `NdnNup`, `NdnNdn` | spin-resolved density-density  | $n_{i\sigma} n_{j\sigma'}$                                       | 2            |
| `SdotS`       | Heisenberg interaction                   | $\mathbf{S}_i \cdot \mathbf{S}_j$                                             | 2            |
| `SzSz`        | Ising interaction                        | $S^z_i S^z_j$                                                                 | 2            |
| `Exchange`    | spin exchange interaction                | $\frac{1}{2}(S^+_i S^-_j + S^-_i S^+_j)$                                      | 2            |
| `ExchangeAsym`| antisymmetric exchange                   | $\frac{1}{2}(S^+_i S^-_j - S^-_i S^+_j)$                                      | 2            |
| `Sz`          | local magnetic moment ($z$)              | $S^z_i$                                                                       | 1            |
| `Sx`, `Sy`    | local magnetic moment ($x$, $y$)         | $S^x_i$, $S^y_i$                                                             | 1            |
| `S+`, `S-`    | spin raising / lowering                  | $S^+_i$, $S^-_i$                                                             | 1            |
| `TotalN`      | total electron number                    | $\sum_i n_i$                                                                  | 0            |
| `TotalNup`, `TotalNdn` | total spin-resolved number      | $\sum_i n_{i\uparrow}$, $\sum_i n_{i\downarrow}$                             | 0            |
| `TotalSz`     | total magnetization                      | $\sum_i S^z_i$                                                                | 0            |
| `Id`          | identity                                 | $\mathbb{1}$                                                                  | 0            |

For a full description of all operator types, see the [operator types](../operators/operator_types.md) page.

## Iteration

An Electron block can be iterated over, where at each iteration a [ProductState](../states/product_state.md) representing the corresponding basis state is returned.

=== "C++"	
	```c++
    auto block = Electron(4, 2, 2);
	for (auto pstate : block) {
      Log("{} {}", to_string(pstate), index(block, pstate));
	}
	```
	
=== "Julia"
	```julia
	block = Electron(4, 2, 2)
	for pstate in block
		@show pstate, index(block, pstate) 
	end
	```

## Methods

#### index

Returns the index of a given [ProductState](../states/product_state.md) in the basis of the tJ block.

=== "C++"	
	```c++
	int64_t index(tJ const &block, ProductState const &pstate);
	```
	
=== "Julia"
	```julia
	index(block::tJ, pstate::ProductState)::Int64
	```
	
!!! warning "1-indexing"
	In the C++ version, the index count starts from "0" whereas in Julia the index count starts from "1".

#### nsites

Returns the number of sites of the block.

=== "C++"	
	```c++
	int64_t nsites(tJ const &block);
	```
	
=== "Julia"
	```julia
	nsites(block::tJ)::Int64
	```

#### size
Returns the size of the block, i.e. its dimension.

=== "C++"	
	```c++
	int64_t size(tJ const &block) const;
	```
	
=== "Julia"
	```julia
	size(block::tJ)::Int64
	```

#### dim
Returns the dimension of the block, same as "size" for non-distributed blocks.

=== "C++"	
	```c++
	int64_t dim(tJ const &block) const;
	```
	
=== "Julia"
	```julia
	dim(block::tJ)::Int64
	```
		
#### isreal
Returns whether the block can be used with real arithmetic. 
Complex arithmetic is needed when a
[Representation](../symmetries/representation.md) is genuinely complex.

=== "C++"	
	```c++
    bool isreal(tJ const &block);
	```

=== "Julia"
	```julia
    isreal(block::tJ)::Bool
	```

## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:Electron"
	```

=== "Julia"
	```julia
	--8<-- "examples/usage_examples/main.jl:Electron"
	```
