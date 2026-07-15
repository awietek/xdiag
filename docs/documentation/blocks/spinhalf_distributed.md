---
title: SpinhalfDistributed
---

A block in a spin $S=1/2$ Hilbert space with distributed computing capabilities.

**Sources:** [spinhalf_distributed.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/blocks/distributed/spinhalf_distributed.hpp) · [spinhalf_distributed.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/blocks/distributed/spinhalf_distributed.cpp)

## Constructors

=== "C++"	
	```c++
	SpinhalfDistributed(int64_t nsites, int64_t nup);
	```
	
| Name    | Description                                                                          | Default |
|:--------|:-------------------------------------------------------------------------------------|---------|
| nsites  | number of sites (integer)                                                            |         |
| nup     | number of "up" spin setting spin (integer)                                           |         |

## Local configurations and operators

A SpinhalfDistributed block describes the same local Hilbert space as the shared-memory [Spinhalf](spinhalf.md) block, only the basis states are distributed across MPI processes. The [local configuration encoding](spinhalf.md#local-configurations) (`0` = ↓, `1` = ↑) is therefore identical to the [Spinhalf](spinhalf.md) block.

Unlike the shared-memory block, the distributed block only supports the operators that have a dedicated distributed kernel. These are:

`Sz`, `S+`, `S-`, `SzSz`, `Exchange`, `ExchangeAsym`, and `Id`.

Their definitions are given in the [Spinhalf operators](spinhalf.md#operators) table.

## Iteration

An SpinhalfDistributed block can be iterated over, where at each iteration a [ProductState](../states/product_state.md) representing the corresponding basis state is returned.

=== "C++"	
	```c++
    auto block = SpinhalfDistributed(4, 2);
	for (auto pstate : block) {
	  Log("{} {}", to_string(pstate), block.index(pstate));
	}
	```
	
## Methods

#### index

Returns the index of a given [ProductState](../states/product_state.md) in the basis of the SpinhalfDistributed block.

=== "C++"	
	```c++
	int64_t SpinhalfDistributed::index(ProductState const &pstate);
	```

#### nsites

Returns the number of sites of the block.

=== "C++"	
	```c++
	int64_t nsites(SpinhalfDistributed const &block);
	```
	
#### size
Returns the size of the block on a local process.

=== "C++"	
	```c++
	int64_t size(SpinhalfDistributed const &block) const;
	```

#### dim
Returns the dimension of the block, i.e. the sum of all sizes across all processes. 

=== "C++"	
	```c++
	int64_t dim(SpinhalfDistributed const &block) const;
	```
		
#### isreal
Returns whether the block can be used with real arithmetic. 
Complex arithmetic is needed when a
[Representation](../symmetries/representation.md) is genuinely complex.

=== "C++"	
	```c++
    bool isreal(SpinhalfDistributed const &block);
	```
