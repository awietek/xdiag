---
title: Permutation
---

# Permutation

## Constructor

Creates an Permutation out of an array of integers, e.g. `[0, 2, 1, 3]`. If the input array is of size `N` then every number between `0` and `N-1` must occur exactly once, otherwise the Permutation is invalid.

=== "Julia"
	```julia
    Permutation(array::Vector{Int64})
	```

=== "C++"	
	```c++
	Permutation(std::vector<int64_t> const &array);
	```

## Methods

### inverse

Computes the inverse permutation.

=== "Julia"
	```julia
    inverse(perm::Permutation)
	```

=== "C++"	
	```c++
	// As a member function
	Permutation inverse() const;
	
	// As a non-member function
	Permutation inverse(Permutation const &p);
	```
---

### "*" operator
Concatenates two permutations by overloading the `*` operator.

=== "Julia"
	```julia
    Base.:*(p1::Permutation, p2::Permutation)
	```

=== "C++"	
	```c++
	Permutation operator*(Permutation const &p1, Permutation const &p2);
	```
---

### size
Returns the size of the permutation, i.e. the number of indices being permuted.

=== "Julia"
	```julia
    size(perm::Permutation)
	```

=== "C++"	
	```c++
	// As a member function
    int64_t size() const;
	
	// As a non-member function
    int64_t size(Permutation const &p);
	```

## Source

[permutation.hpp](https://github.com/awietek/xdiag/blob/master/xdiag/symmetries/permutation.hpp)
