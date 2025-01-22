---
title: Representation
---

A (1D) irreducible representation of a finite group. Upon creation, the group homomorphism properties are verified.

**Source** [representation.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/symmetries/representation.hpp)


## Constructors

### Trivial representation

Creates the trivial representation (all characters equal to 1) of a [PermutationGroup](permutation_group.md)

=== "C++"	
	```c++
	Representation(PermutationGroup const &group);
	```
	
=== "Julia"
	```julia
	Representation(group::PermutationGroup)
	```


	
### With characters

Creates a 1D representation of a [PermutationGroup](permutation_group.md) with given real or complex characters.

=== "C++"	
	```c++
	template <typename T>
	Representation(PermutationGroup const &group, std::vector<T> const &characters);
	template <typename T>
	Representation(PermutationGroup const &group, arma::Col<T> const &characters);
	template <typename T>
	Representation(PermutationGroup const &group, T *characters, int64_t n_characters);
	```
	
| Name  | Description                         |
|:-------------|:---------------------------------------------------------------|
| group        | [PermutationGroup](permutation_group.md) of the Representation |
| characters   | characters of the representation                               |
| n_characters | length of the array of characters                              |

The template parameter `T` in C++ can either be `double` or `complex`.

---

## Methods

#### size
Returns the size of the Representation, i.e. the number of group elements represented.

=== "C++"	
	```c++
	int64_t size(Representation const &irrep);
	```
	
=== "Julia"
	```julia
	size(irrep::Representation)
	```

---
#### isreal 
Returns the whether or not the Representation is real, I.E. the characters are real numbers and do not have an imaginary part.

=== "C++"	
	```c++
	bool isreal(Representation const &irrep) const;
	```
	
=== "Julia"
	```julia
	isreal(irrep::Representation)
	```
---
#### * operator

Multiplies two Representations by overloading the `*` operator.

=== "C++"	
	```c++
	Representation operator*(Representation const &r1, Representation const &r2);
	```

=== "Julia"
	```julia
	Base.:*(r1::Representation, r2::Representation)
	```

	
---

## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:Representation"
	```

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:Representation"
	```
