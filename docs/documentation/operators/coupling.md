---
title: Coupling
---

Describes the coupling of an operator. A coupling can either be a **string** or a [Scalar](scalar.md), which is either a real or complex double precision floating point number. 

**Source** [coupling.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/operators/coupling.hpp)

## Constructors

=== "C++"	
	```c++
	Coupling() = default;
	Coupling(std::string value);
	Coupling(double value);
	Coupling(complex value);
	Coupling(Scalar value););
	```

## Methods

!!! method "isscalar"

	Returns whether or not the coupling is a [Scalar](scalar.md), i.e. a real or complex number.

	=== "C++"	
		```c++
	    bool isscalar(Coupling const& c)
		```
		
!!! method "isstring"

	Returns whether or not the coupling is a string.

	=== "C++"	
		```c++
	    bool isstring(Coupling const& c)
		```
		
!!! method "scalar"

	Returns the [Scalar](scalar.md) if the Coupling holds a scalar. Otherwise throws an `xdiag::Error`

	=== "C++"	
		```c++
	    Scalar scalar(Coupling const& c)
		```
!!! method "string"

	Returns the string if the Coupling holds a string. Otherwise throws an `xdiag::Error`

	=== "C++"	
		```c++
	    std::string string(Coupling const& c)
		```
