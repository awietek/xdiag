---
title: PermutationGroup
---

A group of permutations. Group axioms are verified during construction.

**Source** [permutation_group.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/symmetries/permutation_group.hpp)

## Constructor

### Constructor from Permutations

Creates an PermutationGroup out of a vector of [Permutation](permutation.md) objects.


=== "C++"	
	```c++
	PermutationGroup(std::vector<Permutation> const &permutations);
	```
	
=== "Julia"
	```julia
	PermutationGroup(permutations::Vector{Permutation})
	```


### Constructor from Matrix

Creates a PermutationGroup out of a matrix whose rows specify the individual permutations.

=== "C++"	
	```c++
    PermutationGroup(arma::Mat<int64_t> const &matrix);
	PermutationGroup(int64_t *ptr, int64_t n_permutations, int64_t nsites);
	```
---

## Methods


#### nsites

Returns the number of sites on which the permutations of the group acts.

=== "C++"	
	```c++
	int64_t nsites(PermutationGroup const &group);
	```
=== "Julia"
	```julia
	nsites(group::PermutationGroup)
	```

---

#### size
Returns the size of the permutation group, i.e. the number permutations.

=== "C++"	
	```c++
	int64_t size(PermutationGroup const &group);
	```
=== "Julia"
	```julia
	size(group::PermutationGroup)
	```
---

## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:PermutationGroup"
	```
=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:PermutationGroup"
	```
