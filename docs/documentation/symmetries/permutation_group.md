---
title: PermutationGroup
---

A group of permutations. Group axioms are verified during construction.

**Sources**<br>
[permutation_group.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/symmetries/permutation_group.hpp)<br>
[permutation_group.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/symmetries/permutation_group.cpp)<br>
[permutation_group.jl](https://github.com/awietek/XDiag.jl/blob/main/src/symmetries/permutation_group.jl)


---

## Constructor

### From Permutations

Creates an PermutationGroup out of a vector of [Permutation](permutation.md) objects.


=== "C++"	
	```c++
	PermutationGroup(std::vector<Permutation> const &permutations);
	```
	
=== "Julia"
	```julia
	PermutationGroup(permutations::Vector{Permutation})
	```


### From matrix

Creates a PermutationGroup out of a matrix whose rows specify the individual permutations. If a raw pointer is handed, the matrix is assumed to be in column-major form.

=== "C++"	
	```c++
    PermutationGroup(arma::Mat<int64_t> const &matrix);
	PermutationGroup(int64_t *ptr, int64_t n_permutations, int64_t nsites);
	```
=== "Julia"
	```julia
	PermutationGroup(matrix::Matrix{Int64})
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
	nsites(group::PermutationGroup)::Int64
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
	size(group::PermutationGroup)::Int64
	```
---

#### to_string (operator<<)

Converts the PermutationGroup to a readable string representation.
	
=== "C++"	
	```c++
	std::string to_string(PermutationGroup const &group);
	std::ostream &operator<<(std::ostream &out, PermutationGroup const &group);
	```

=== "Julia"
	```julia
    to_string(group::PermutationGroup)::String
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
