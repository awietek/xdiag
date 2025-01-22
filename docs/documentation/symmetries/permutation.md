---
title: Permutation
---

Permutations of indices or lattice sites. Basic building block of a [PermutationGroup](permutation_group.md). Permutations can be multiplied, inverted and raised to a power.

**Source** [permutation.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/symmetries/permutation.hpp)

## Constructors

	
### Constructors from an array

Creates an Permutation out of an array of integers, e.g. `{0, 2, 1, 3}`. If the input array is of size `N` then every number between `0` and `N-1` must occur exactly once, otherwise the Permutation is invalid.

!!! warning "1-indexing in Julia / 0-indexing in C++"

	To enumerate the sites of a Permutation, we start counting at 1 in Julia and 0 in C++.
	

=== "C++"	
	```c++
	Permutation(std::initializer_list<int64_t> list);
	Permutation(std::vector<int32_t> const &array);
	Permutation(std::vector<int64_t> const &array);
	Permutation(arma::Col<int64_t> const &array);
	Permutation(int64_t *array, int64_t size);
	```
=== "Julia"
	```julia
 	Permutation(array::Vector{Int64})
 	```
	
| Name  | Description                         |
|:------|:------------------------------------|
| array | array of integers, e.g. {0,2,1,3}   |
| list  | initializer list of the permutation |
| ptr   | pointer to memory as an array       |
| size  | size of the array                   |



### Constructor for Identity

Constructs an identity permutation of a given size, e.g. `{0, 1, 2, 3}`.

=== "C++"	
	```c++
	Permutation(int64_t size);
	```
=== "Julia"
	```julia
 	Permutation(size::Integer)
 	```
	
| Name | Description                      |
|:-----|:---------------------------------|
| size | size of the identity permutation |

---

## Methods


#### inverse

Computes the inverse permutation.
	
=== "C++"	
	```c++
	Permutation inverse(Permutation const &p);
	```
	
=== "Julia"
	```julia
	inverse(perm::Permutation)
	```
---

#### * operator

Concatenates two permutations by overloading the `*` operator.

=== "C++"	
	```c++
	Permutation operator*(Permutation const &p1, Permutation const &p2);
	```
	
=== "Julia"
	```julia
	Base.:*(p1::Permutation, p2::Permutation)
	```
	
---
#### ^ operator, pow

Raises a permutation to an integer power.


=== "C++"	
	```c++
	Permutation pow(Permutation const &p, int64_t power);
	```
=== "Julia"
	```julia
	Base.:^(p::Permutation, power::Integer)
	```
---
#### size

Returns the size of a Permutation.

=== "C++"	
	```c++
	int64_t size(Permutation const &p);
	```
=== "Julia"
	```julia
	size(p::Permutation)
	```

## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:Permutation"
	```

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:Permutation"
	```
