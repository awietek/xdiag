---
title: Electron
---

# Electron

## Constructors

Creates a block object for Electron (spinful fermion) type Hilbert spaces. 


???+ method "Parameters"

    | Name    | Description                                                                                                     |
    |:--------|:----------------------------------------------------------------------------------------------------------------|
    | n_sites | number of sites (integer)                                                                                       |
    | n_up    | number of "up" electrons (integer)                                                                              |
    | n_dn    | number of "dn" electrons (integer)                                                                              |
    | group   | [PermutationGroup](../symmetries/permutation_group.md) defining the permutation symmetries                      |
    | irrep   | [Representation](../symmetries/representation.md) defining the irreducible representation of the symmetry group |

=== "Julia"
	```julia
	Electron(n_sites::Integer)
	Electron(n_sites::Integer, n_up::Integer, n_dn::Integer)
	Electron(n_sites::Integer, group::PermutationGroup, irrep::Representation)
	Electron(n_sites::Integer, n_up::Integer, n_dn::Integer, 
	         group::PermutationGroup, irrep::Representation)
	```

=== "C++"	
	```c++
    Electron(int64_t n_sites);
    Electron(int64_t n_sites, int64_t n_up, int64_t n_dn);
    Electron(int64_t n_sites, PermutationGroup permutation_group,
             Representation irrep);
    Electron(int64_t n_sites, int64_t n_up, int64_t n_dn, 
	         PermutationGroup group, Representation irrep);
	```


## Methods


??? method "n_sites"

	Returns the number of sites of the block.

	=== "Julia"
		```julia
		n_sites(block::Electron)
		```

	=== "C++"	
		```c++
		int64_t n_sites() const
		```

??? method "size"
	Returns the size of the block, i.e. its dimension.

	=== "Julia"
		```julia
		size(block::Electron)
		```

	=== "C++"	
		```c++
		int64_t size() const;
		```


## Usage Example

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:Electron"
	```

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:Electron"
	```

**Source** [spinhalf.hpp](https://github.com/awietek/xdiag/blob/master/xdiag/blocks/electron/electron.hpp)
