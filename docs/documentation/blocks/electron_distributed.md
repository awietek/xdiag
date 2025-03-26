---
title: ElectronDistributed
---

A block in an electron type Hilbert space, i.e. fermions with $\uparrow, \downarrow$ spin with distributed computing capabilities. 

**Sources**<br>
[electron_distributed.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/blocks/electron_distributed.hpp)<br>
[electron_distributed.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/blocks/electron_distributed.cpp)

---

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
| backend | backend used for coding the basis states                                             | `auto`  |

The parameter `backend` chooses how the block is coded internally. By using the default parameter `auto` the backend is chosen automatically. Alternatives are `32bit`, `64bit`.

---

## Iteration

An ElectronDistributed block can be iterated over, where at each iteration a [ProductState](../states/product_state.md) representing the corresponding basis state is returned.

=== "C++"	
	```c++
    auto block = ElectronDistributed(4, 2, 1);
	for (auto pstate : block) {
		Log("{} {}", to_string(pstate), index(block, pstate));
	}
	```

---

## Methods

#### index

Returns the index of a given [ProductState](../states/product_state.md) in the basis of the ElectronDistributed block.

=== "C++"	
	```c++
	int64_t index(ElectronDistributed const &block, ProductState const &pstate);
	```

---

#### nsites

Returns the number of sites of the block.

=== "C++"	
	```c++
	int64_t nsites(ElectronDistributed const &block);
	```
---

#### size
Returns the size of the block on a local process.

=== "C++"	
	```c++
	int64_t size(ElectronDistributed const &block) const;
	```

---

#### dim
Returns the dimension of the block, i.e. the sum of all sizes across all processes. 

=== "C++"	
	```c++
	int64_t dim(ElectronDistributed const &block) const;
	```
---
		
#### isreal
Returns whether the block can be used with real arithmetic. 
Complex arithmetic is needed when a
[Representation](../symmetries/representation.md) is genuinely complex.

=== "C++"	
	```c++
    bool isreal(ElectronDistributed const &block);
	```
