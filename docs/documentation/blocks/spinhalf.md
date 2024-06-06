---
title: Spinhalf
---

# Spinhalf

## Constructors

Creates a block object for spin $S=1/2$ type Hilbert spaces. 


???+ method "Parameters"

    | Name    | Description                                                                                                     |
    |:--------|:----------------------------------------------------------------------------------------------------------------|
    | n_sites | number of sites (integer)                                                                                       |
    | n_up    | number of "up" spin setting spin (integer)                                                                      |
    | group   | [PermutationGroup](../symmetries/permutation_group.md) defining the permutation symmetries                      |
    | irrep   | [Representation](../symmetries/representation.md) defining the irreducible representation of the symmetry group |

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


## Methods


??? method "n_sites"

	Returns the number of sites of the block.

	=== "Julia"
		```julia
		n_sites(block::Spinhalf)
		```

	=== "C++"	
		```c++
		int64_t n_sites() const
		```

??? method "size"
	Returns the size of the block, i.e. its dimension.

	=== "Julia"
		```julia
		size(block::Spinhalf)
		```

	=== "C++"	
		```c++
		int64_t size() const;
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

**Source** [spinhalf.hpp](https://github.com/awietek/xdiag/blob/master/xdiag/blocks/spinhalf/spinhalf.hpp)
