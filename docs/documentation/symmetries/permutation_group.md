---
title: Permutation
---

# PermutationGroup

## Constructor

Creates an PermutationGroup out of a vector of Permutations.

=== "Julia"
	```julia
	PermutationGroup(permutations::Vector{Permutation})
	```

=== "C++"	
	```c++
	PermutationGroup(std::vector<Permutation> const &permutations);
	```

## Methods


??? method "n_sites"

	Returns the number of sites on which the permutations of the group acts.

	=== "Julia"
		```julia
		n_sites(group::PermutationGroup)
		```

	=== "C++"	
		```c++
		int64_t n_sites() const
		```

??? method "size"
	Returns the size of the permutation group, i.e. the number permutations

	=== "Julia"
		```julia
		size(group::PermutationGroup)
		```

	=== "C++"	
		```c++
		int64_t size() const;
		```

??? method "inverse"

	Given an index of a permutation, it returns the index of the inverse permutation.

	=== "Julia"
		```julia
		inverse(group::PermutationGroup, idx::Integer)
		```

	=== "C++"	
		```c++
		// As a member function
        int64_t inverse(int64_t sym) const;
		```


## Usage Example

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:PermutationGroup"
	```

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:PermutationGroup"
	```

**Source** [permutation_group.hpp](https://github.com/awietek/xdiag/blob/master/xdiag/symmetries/permutation_group.hpp)
