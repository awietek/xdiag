---
title: PermutationGroup
---

A group of permutations. Group axioms are verified during construction.

**Source** [permutation_group.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/symmetries/permutation_group.hpp)

## Constructor

### Constructor from Permutations

Creates an PermutationGroup out of a vector of [Permutation](permutation.md) objects.

=== "Julia"
	```julia
	PermutationGroup(permutations::Vector{Permutation})
	```

=== "C++"	
	```c++
	PermutationGroup(std::vector<Permutation> const &permutations);
	```

### Constructor from Matrix

Creates a PermutationGroup out of a matrix whose rows specify the individual permutations.

=== "C++"	
	```c++
    PermutationGroup(arma::Mat<int64_t> const &matrix);
	PermutationGroup(int64_t *ptr, int64_t n_permutations, int64_t n_sites);
	```
---

## Methods


#### n_sites

Returns the number of sites on which the permutations of the group acts.

=== "Julia"
	```julia
	n_sites(group::PermutationGroup)
	```

=== "C++"	
	```c++
	int64_t n_sites(PermutationGroup const &group);
	```
---

#### size
Returns the size of the permutation group, i.e. the number permutations.
	
=== "Julia"
	```julia
	size(group::PermutationGroup)
	```

=== "C++"	
	```c++
	int64_t size(PermutationGroup const &group);
	```
---

## Usage Example

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:PermutationGroup"
	```

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:PermutationGroup"
	```
