---
title: tJ
---

A block in a  $t-J$ type Hilbert space, i.e. fermions with $\uparrow, \downarrow$ spin excluding doubly occupied sites. 

**Sources**<br>
[tj.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/blocks/tj.hpp)<br>
[tj.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/blocks/tj.cpp)<br>
[tj.jl](https://github.com/awietek/XDiag.jl/blob/main/src/blocks/tj.jl)
---

## Constructors

=== "C++"	
	```c++
	tJ(int64_t nsites, int64_t nup, int64_t ndn, std::string backend = "auto");
	tJ(int64_t nsites, int64_t nup, int64_t ndn, Representation const &irrep, std::string backend = "auto");
	```

=== "Julia"
	```julia
	tJ(nsites::Int64, nup::Int64, ndn::Int64, backend::String="auto")
	tJ(nsites::Int64, nup::Int64, ndn::Int64, irrep::Representation, backend::String="auto")
	```


| Name    | Description                                                                          | Default |
|:--------|:-------------------------------------------------------------------------------------|---------|
| nsites  | number of sites (integer)                                                            |         |
| nup     | number of "up" electrons (integer)                                                   |         |
| ndn     | number of "dn" electrons (integer)                                                   |         |
| irrep   | Irreducible [Representation](../symmetries/representation.md)  of the symmetry group |         |
| backend | backend used for coding the basis states                                             | `auto`  |

The parameter `backend` chooses how the block is coded internally. By using the default parameter `auto` the backend is chosen automatically. Alternatives are `32bit`, `64bit`.

---

## Iteration

An tJ block can be iterated over, where at each iteration a [ProductState](../states/product_state.md) representing the corresponding basis state is returned.

=== "C++"	
	```c++
    auto block = tJ(4, 2, 1);
	for (auto pstate : block) {
		Log("{} {}", to_string(pstate), index(block, pstate));
	}
	```
=== "Julia"
	```julia
	block = tJ(4, 2, 1)
	for pstate in block
		@show pstate, index(block, pstate) 
	end
	```
---

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

---

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
---

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

---

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

---
		
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
    isreal(block::tJ)::Int64
	```
---

## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:tJ"
	```

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:tJ"
	```

