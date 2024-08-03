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


## Methods


!!! method "n_sites"

	Returns the number of sites of the block.

	=== "Julia"
		```julia
		n_sites(block::Spinhalf)
		```

	=== "C++"	
		```c++
		int64_t n_sites() const
		```

!!! method "dim"
	Returns the dimension of the block.

	=== "Julia"
		```julia
		dim(block::Spinhalf)
		```

	=== "C++"	
		```c++
		int64_t dim() const;
		```


!!! method "size"
	Returns the size of the block locally. Same as "dim" for non-distributed Blocks but different for distributed blocks.

	=== "Julia"
		```julia
		size(block::Spinhalf)
		```

	=== "C++"	
		```c++
		int64_t size() const;
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

