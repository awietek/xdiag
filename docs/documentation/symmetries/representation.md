---
title: Representation
---

## Constructors

Creates a Representation from a vector of complex numbers

=== "Julia"
	```julia
	Representation(characters::Vector{<:Number})
	```

=== "C++"	
	```c++
	Representation(std::vector<complex> const &characters);
	```

## Methods

??? method "size"
	Returns the size of the Representation, i.e. the number of characters.

	=== "Julia"
		```julia
		size(irrep::Representation)
		```

	=== "C++"	
		```c++
        int64_t size() const;
		```


??? method "isreal"
	Returns the whether or not the Representation is real.

	=== "Julia"
		```julia
		isreal(irrep::Representation; precision::Float64=1e-12)
		```

	=== "C++"	
		```c++
        bool isreal(double precision = 1e-12) const;
		```

??? method ""*" operator"

	Multiplies two Representations by overloading the `*` operator.

	=== "Julia"
		```julia
		Base.:*(p1::Representation, p2::Representation)
		```

	=== "C++"	
		```c++
		Representation operator*(Representation const &p1, Representation const &p2);
		```

## Usage Example

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:Representation"
	```

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:Representation"
	```

**Source** [representation.hpp](https://github.com/awietek/xdiag/blob/master/xdiag/symmetries/representation.hpp)
