---
title: Spinhalf
---

A block in a spin $S=1/2$  Hilbert space. 

**Source** [spinhalf.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/blocks/spinhalf.hpp)


## Constructors

=== "C++"	
	```c++
	Spinhalf(int64_t nsites, std::string backend = "auto");
	Spinhalf(int64_t nsites, int64_t nup, std::string backend = "auto");
	Spinhalf(int64_t nsites, Representation const &irrep, std::string backend = "auto");
	Spinhalf(int64_t nsites, int64_t nup, Representation const &irrep, std::string backend = "auto");

	```
=== "Julia"
	```julia
	Spinhalf(nsites::Integer, backend::AbstractString)
	Spinhalf(nsites::Integer, nup::Integer, backend::AbstractString)
	Spinhalf(nsites::Integer, irrep::Representation, backend::AbstractString)
	Spinhalf(nsites::Integer, nup::Integer, irrep::Representation, backend::AbstractString)
	```
	
| Name    | Description                                                                          | Default |
|:--------|:-------------------------------------------------------------------------------------|---------|
| nsites  | number of sites (integer)                                                            |         |
| nup     | number of "up" spin setting spin (integer)                                           |         |
| irrep   | Irreducible [Representation](../symmetries/representation.md)  of the symmetry group |         |
| backend | backend used for coding the basis states                                             | `auto`  |
	
	
The parameter `backend` chooses how the block is coded internally. By using the default parameter `auto` the backend is chosen automatically. Alternatives are `32bit`, `64bit`, `1sublattice`, `2sublattice`, `3sublattice`, `4sublattice`, and `5sublattice`. The backends `xsublattice` implement the sublattice coding algorithm described in [Wietek, Läuchli, Phys. Rev. E 98, 033309 (2018)](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.98.033309). The sublattice coding algorithms impose certain constraints on the symmetries used, as described in the reference. 

## Iteration

An Spinhalf block can be iterated over, where at each iteration a [ProductState](../states/product_state.md) representing the corresponding basis state is returned.

=== "C++"	
	```c++
    auto block = Spinhalf(4, 2);
	for (auto pstate : block) {
	  Log("{} {}", to_string(pstate), block.index(pstate));
	}
	```
=== "Julia"
	```julia
	block = Spinhalf(4, 2)
	for pstate in block
		@show pstate, index(block, pstate) 
	end
	```


	
	
---

## Methods

#### index

Returns the index of a given [ProductState](../states/product_state.md) in the basis of the Spinhalf block.

=== "C++"	
	```c++
	int64_t index(Spinhalf const &block, ProductState const &pstate);
	```
	
=== "Julia"
	```julia
	index(block::Spinhalf, pstate::ProductState)
	```
	
!!! warning "1-indexing"
	In the C++ version, the index count starts from "0" whereas in Julia the index count starts from "1".

---

#### nsites

Returns the number of sites of the block.

=== "C++"	
	```c++
	int64_t nsites(Spinhalf const &block);
	```
	
=== "Julia"
	```julia
	nsites(block::Spinhalf)
	```
---

#### size
Returns the size of the block, i.e. its dimension.

=== "C++"	
	```c++
	int64_t size(Spinhalf const &block) const;
	```
	
=== "Julia"
	```julia
	size(block::Spinhalf)
	```

---

#### dim
Returns the dimension of the block, same as "size" for non-distributed blocks.

=== "C++"	
	```c++
	int64_t dim(Spinhalf const &block) const;
	```
	
=== "Julia"
	```julia
	dim(block::Spinhalf)
	```

---
		
#### isreal
Returns whether the block can be used with real arithmetic. 
Complex arithmetic is needed when a
[Representation](../symmetries/representation.md) is genuinely complex.

=== "C++"	
	```c++
    bool isreal(Spinhalf const &block);
	```

=== "Julia"
	```julia
    isreal(block::Spinhalf)
	```
---

## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:Spinhalf"
	```

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:Spinhalf"
	```
