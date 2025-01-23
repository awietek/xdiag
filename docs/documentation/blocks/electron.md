---
title: Electron
---

Representation of a block in a Electron (fermions with $\uparrow, \downarrow$ spin) Hilbert space. 

**Source** [electron.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/blocks/electron.hpp)

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
	Electron(nsites::Integer, backend::String)
	Electron(nsites::Integer, nup::Integer, ndn::Integer, backend::String)
	Electron(nsites::Integer, irrep::Representation, backend::String)
	Electron(nsites::Integer, nup::Integer, ndn::Integer, irrep::Representation, backend::String)
	```

| Name    | Description                                                                          |        |
|:--------|:-------------------------------------------------------------------------------------|--------|
| nsites  | number of sites (integer)                                                            |        |
| nup     | number of "up" electrons (integer)                                                   |        |
| ndn     | number of "dn" electrons (integer)                                                   |        |
| irrep   | Irreducible [Representation](../symmetries/representation.md)  of the symmetry group |        |
| backend | backend used for coding the basis states                                             | `auto` |

The parameter `backend` chooses how the block is coded internally. By using the default parameter `auto` the backend is chosen automatically. Alternatives are `32bit`, `64bit`

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
	index(block::tJ, pstate::ProductState)
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
	nsites(block::tJ)
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
	size(block::tJ)
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
	dim(block::tJ)
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
    isreal(block::tJ)
	```
---

## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:Electron"
	```

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:Electron"
	```
