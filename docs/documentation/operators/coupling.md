---
title: Coupling
---

Describes the coupling of a local operator. A coupling can either be a string, a real/complex number or even a real/complex matrix. It allows for converting to real/complex numbers or matrices as well as strings, whenever this conversion is sensible. 

**Source** [coupling.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/operators/coupling.hpp)

## Constructors

=== "Julia"
	```julia
	Coupling(name::String)
	Coupling(val::Float64)
	Coupling(val::ComplexF64)
	Coupling(mat::Matrix{Float64})
   	Coupling(mat::Matrix{ComplexF64})
	```
=== "C++"	
	```c++
	Coupling(std::string value);
	Coupling(double value);
	Coupling(complex value);
	Coupling(arma::mat const &value);
	Coupling(arma::cx_mat const &value);
	```

## Methods

!!! method "type"

	Returns the type of the Coupling, i.e. a string which either reads "string", "double", "complex", "mat", or "cx_mat"

	=== "Julia"
		```julia
		type(cpl::Coupling)
		```

	=== "C++"	
		```c++
		std::string type() const;
		```
		
!!! method "isreal"

	Returns whether or not the coupling is real. Throws an error if the coupling is given as a string, since then it cannot be determined whether the operator is real.

	=== "Julia"
		```julia
		isreal(cpl::Coupling)
		```

	=== "C++"	
		```c++
		bool isreal() const;
		```	
		
!!! method "ismatrix"

	Returns whether or not the coupling is defined as a matrix. Throws an error if the coupling is given as a string, since then it cannot be determined whether the operator is real.

	=== "Julia"
		```julia
		ismatrix(cpl::Coupling)
		```

	=== "C++"	
		```c++
		bool ismatrix() const;
		```	
		
				
!!! method "isexplicit"

	Returns **false** if the coupling is defined as a string, otherwise **true**

	=== "Julia"
		```julia
		isexplicit(cpl::Coupling)
		```

	=== "C++"	
		```c++
		bool isexplicit() const;
		```	

## Conversions

A Coupling can be converted to the values it represents, so a string, real/complex number or a real/complex matrix. Initially real values can be cast to complex.

=== "Julia"
	```julia
	convert(::Type{String}, cpl::Coupling)
	convert(::Type{Float64}, cpl::Coupling)
	convert(::Type{ComplexF64}, cpl::Coupling)
	convert(::Type{Matrix{Float64}}, cpl::Coupling)
	convert(::Type{Matrix{ComplexF64}}, cpl::Coupling)
	```
=== "C++"	
	```c++
	template <typename coeff_t> coeff_t as() const;
	```

## Usage Example

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:coupling"
	```

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:coupling"
	```


