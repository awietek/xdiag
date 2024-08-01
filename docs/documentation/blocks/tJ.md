---
title: tJ
---

Representation of a block in a  $t-J$ type Hilbert space. 

**Source** [tj.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/blocks/tj.hpp)

## Constructors

=== "Julia"
	```julia
	tJ(n_sites::Integer, n_up::Integer, n_dn::Integer)
	tJ(n_sites::Integer, n_up::Integer, n_dn::Integer, 
	   group::PermutationGroup, irrep::Representation)
	```

=== "C++"	
	```c++
    tJ(int64_t n_sites, int64_t n_up, int64_t n_dn);
    tJ(int64_t n_sites, int64_t n_up, int64_t n_dn, 
	   PermutationGroup group, Representation irrep);
	```


| Name    | Description                                                                                |   |
|:--------|:-------------------------------------------------------------------------------------------|---|
| n_sites | number of sites (integer)                                                                  |   |
| n_up    | number of "up" electrons (integer)                                                         |   |
| n_dn    | number of "dn" electrons (integer)                                                         |   |
| group   | [PermutationGroup](../symmetries/permutation_group.md) defining the permutation symmetries |   |
| irrep   | Irreducible [Representation](../symmetries/representation.md)  of the symmetry group       |   |


## Methods


!!! method "n_sites"

	Returns the number of sites of the block.

	=== "Julia"
		```julia
		n_sites(block::tJ)
		```

	=== "C++"	
		```c++
		int64_t n_sites() const
		```

!!! method "size"
	Returns the size of the block, i.e. its dimension.

	=== "Julia"
		```julia
		size(block::tJ)
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
	    isreal(block::tJ; precision::Real=1e-12)
		```

	=== "C++"	
		```c++
		int64_t isreal(double precision = 1e-12) const;
		```


## Usage Example

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:tJ"
	```

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:tJ"
	```

