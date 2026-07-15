---
title: ElectronDistributed
---

A block in an electron type Hilbert space, i.e. fermions with $\uparrow, \downarrow$ spin with distributed computing capabilities. 

**Sources:** [electron_distributed.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/blocks/distributed/electron_distributed.hpp) · [electron_distributed.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/blocks/distributed/electron_distributed.cpp)

## Constructors

=== "C++"	
	```c++
	ElectronDistributed(int64_t nsites, int64_t nup, int64_t ndn, std::string backend = "auto");
	```

| Name    | Description                                                                          | Default |
|:--------|:-------------------------------------------------------------------------------------|---------|
| nsites  | number of sites (integer)                                                            |         |
| nup     | number of "up" electrons (integer)                                                   |         |
| ndn     | number of "dn" electrons (integer)                                                   |         |

## Local configurations and operators

An ElectronDistributed block describes the same local Hilbert space as the shared-memory [Electron](electron.md) block, only the basis states are distributed across MPI processes. The [local configuration encoding](electron.md#local-configurations) (`0` = empty, `1` = ↑, `2` = ↓, `3` = ⇅) is therefore identical to the [Electron](electron.md) block.

Unlike the shared-memory block, the distributed block only supports the operators that have a dedicated distributed kernel. These are:

`Cdagup`, `Cup`, `Cdagdn`, `Cdn`, `Hopup`, `Hopdn`, `HopupAsym`, `HopdnAsym`, `Nup`, `Ndn`, `Nupdn`, `NupNup`, `NupNdn`, `NdnNup`, `NdnNdn`, `NtotNtot`, `NupdnNupdn`, `HubbardU`, `SzSz`, `Exchange`, `ExchangeAsym`, and `Id`.

Their definitions are given in the [Electron operators](electron.md#operators) table.

## Iteration

An ElectronDistributed block can be iterated over, where at each iteration a [ProductState](../states/product_state.md) representing the corresponding basis state is returned.

=== "C++"	
	```c++
    auto block = ElectronDistributed(4, 2, 1);
	for (auto pstate : block) {
		Log("{} {}", to_string(pstate), block.index(pstate));
	}
	```

## Methods

#### index

Returns the index of a given [ProductState](../states/product_state.md) in the basis of the ElectronDistributed block.

=== "C++"	
	```c++
	int64_t ElectronDistributed::index(ProductState const &pstate);
	```

#### nsites

Returns the number of sites of the block.

=== "C++"	
	```c++
	int64_t nsites(ElectronDistributed const &block);
	```

#### size
Returns the size of the block on a local process.

=== "C++"	
	```c++
	int64_t size(ElectronDistributed const &block) const;
	```

#### dim
Returns the dimension of the block, i.e. the sum of all sizes across all processes. 

=== "C++"	
	```c++
	int64_t dim(ElectronDistributed const &block) const;
	```
		
#### isreal
Returns whether the block can be used with real arithmetic. 
Complex arithmetic is needed when a
[Representation](../symmetries/representation.md) is genuinely complex.

=== "C++"	
	```c++
    bool isreal(ElectronDistributed const &block);
	```
