---
title: Permutation
---

Permutations of indices or lattice sites

**Source** [permutation.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/symmetries/permutation.hpp)

## Constructors

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

!!! method "inverse"

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

!!! method ""*" operator"

	Concatenates two permutations by overloading the `*` operator.

	=== "Julia"
		```julia
		Base.:*(p1::Permutation, p2::Permutation)
		```

	=== "C++"	
		```c++
		Permutation operator*(Permutation const &p1, Permutation const &p2);
		```

!!! method "size"
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

## Usage Example

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:Permutation"
	```

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:Permutation"
	```

