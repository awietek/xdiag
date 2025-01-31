---
title: tJDistributed
---

A block in a  $t-J$ type Hilbert space, i.e. fermions with $\uparrow, \downarrow$ spin excluding doubly occupied sites with distributed computing capabilities. 

**Sources** [tj_distributed.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/blocks/tj_distributed.hpp), [tj_distributed.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/blocks/tj_distributed.cpp)

---

## Constructors

=== "C++"	
	```c++
	tJDistributed(int64_t nsites, int64_t nup, int64_t ndn, std::string backend = "auto");
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

An tJDistributed block can be iterated over, where at each iteration a [ProductState](../states/product_state.md) representing the corresponding basis state is returned.

=== "C++"	
	```c++
    auto block = tJDistributed(4, 2, 1);
	for (auto pstate : block) {
		Log("{} {}", to_string(pstate), index(block, pstate));
	}
	```

---

## Methods

#### index

Returns the index of a given [ProductState](../states/product_state.md) in the basis of the tJDistributed block.

=== "C++"	
	```c++
	int64_t index(tJDistributed const &block, ProductState const &pstate);
	```

---

#### nsites

Returns the number of sites of the block.

=== "C++"	
	```c++
	int64_t nsites(tJDistributed const &block);
	```

---

#### size
Returns the size of the block on a local process.

=== "C++"	
	```c++
	int64_t size(tJDistributed const &block) const;
	```


---

#### dim
Returns the dimension of the block, i.e. the sum of all sizes across all processes. 

=== "C++"	
	```c++
	int64_t dim(tJDistributed const &block) const;
	```


---
		
#### isreal
Returns whether the block can be used with real arithmetic. 
Complex arithmetic is needed when a
[Representation](../symmetries/representation.md) is genuinely complex.

=== "C++"	
	```c++
    bool isreal(tJDistributed const &block);
	```
