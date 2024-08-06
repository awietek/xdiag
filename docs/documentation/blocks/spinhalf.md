---
title: Spinhalf
---

Representation of a block in a spin $S=1/2$  Hilbert space. 

**Source** [spinhalf.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/blocks/spinhalf.hpp)


## Constructors
=== "Julia"
	```julia
	Spinhalf(n_sites::Integer)
	Spinhalf(n_sites::Integer, n_up::Integer)
	Spinhalf(n_sites::Integer, group::PermutationGroup, irrep::Representation)
	Spinhalf(n_sites::Integer, n_up::Integer, group::PermutationGroup, 
	         irrep::Representation)
	```

=== "C++"	
	```c++
    Spinhalf(int64_t n_sites);
    Spinhalf(int64_t n_sites, int64_t n_up);
    Spinhalf(int64_t n_sites, PermutationGroup permutation_group,
             Representation irrep);
    Spinhalf(int64_t n_sites, int64_t n_up, PermutationGroup group,
             Representation irrep);
	```
	
| Name    | Description                                                                                |   |
|:--------|:-------------------------------------------------------------------------------------------|---|
| n_sites | number of sites (integer)                                                                  |   |
| n_up    | number of "up" spin setting spin (integer)                                                 |   |
| group   | [PermutationGroup](../symmetries/permutation_group.md) defining the permutation symmetries |   |
| irrep   | Irreducible [Representation](../symmetries/representation.md)  of the symmetry group       |   |


## Iteration

An Spinhalf block can be iterated over, where at each iteration a [ProductState](../states/product_state.md) representing the corresponding basis state is returned.

=== "Julia"
	```julia
	block = Spinhalf(4, 2)
	for pstate in block
		@show pstate, index(block, pstate) 
	end
	```

=== "C++"	
	```c++
    auto block = Spinhalf(4, 2);
	for (auto pstate : block) {
		Log("{} {}", to_string(pstate), block.index(pstate));
	}
	```
	
	

## Methods

!!! method "index"

	Returns the index of a given [ProductState](../states/product_state.md) in the basis of the Spinhalf block.

	=== "Julia"
		```julia
		index(block::Spinhalf, pstate::ProductState)
		```

	=== "C++"	
		```c++
		int64_t index(ProductState const &pstate) const;
		```
		
	!!! warning "1-indexing"
		In the C++ version, the index count starts from "0" whereas in Julia the index count starts from "1".



!!! method "n_sites"

	Returns the number of sites of the block.

	=== "Julia"
		```julia
		n_sites(block::Spinhalf)
		```

	=== "C++"	
		```c++
		int64_t n_sites() const;
		```

!!! method "n_up"

	Returns the number of "up" spins.

	=== "Julia"
		```julia
		n_up(block::Spinhalf)
		```

	=== "C++"	
		```c++
		int64_t n_up() const;
		```

!!! method "permutation_group"

	Returns the [PermutationGroup]("../symmetries/permutation_group.md") of the block, if defined.

	=== "Julia"
		```julia
		permutation_group(block::Spinhalf)
		```

	=== "C++"	
		```c++
	    PermutationGroup permutation_group() const;
		```


!!! method "irrep"

	Returns the [Representation]("../symmetries/representation.md") of the block, if defined.

	=== "Julia"
		```julia
	    irrep(block::Spinhalf)
		```

	=== "C++"	
		```c++
	    Representation irrep() const;
		```


!!! method "size"
	Returns the size of the block, i.e. its dimension.

	=== "Julia"
		```julia
		size(block::Spinhalf)
		```

	=== "C++"	
		```c++
		int64_t size() const;
		```

!!! method "dim"
	Returns the dimension of the block, same as "size" for non-distributed blocks.

	=== "Julia"
		```julia
		dim(block::Spinhalf)
		```

	=== "C++"	
		```c++
		int64_tdim() const;
		```
		
!!! method "isreal"
	Returns whether the block can be used with real arithmetic. 
	Complex arithmetic is needed when a
	[Representation](../symmetries/representation.md) is genuinely complex.

	=== "Julia"
		```julia
	    isreal(block::Spinhalf; precision::Real=1e-12)
		```

	=== "C++"	
		```c++
		int64_t isreal(double precision = 1e-12) const;
		```


## Usage Example

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:Spinhalf"
	```

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:Spinhalf"
	```

