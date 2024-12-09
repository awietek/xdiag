---
title: Scalar
---

A variable which can be either a real or complex double precision number. The algebraic operators `+,-,*,\` and comparison operators `==, !=` are defined.

**Source** [coupling.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/operators/scalar.hpp)

## Constructors

=== "C++"	
	```c++
	Scalar() = default;
	explicit Scalar(double value);
	explicit Scalar(complex value);
	```
## Methods

!!! method "isreal"

    Returns `true` if Scalar is a real number, `false` otherwise. 
    === "C++" 
    	```c++
		bool isreal(Scalar const &s)
		```

!!! method "real"

    Returns the real part of the Scalar.
    === "C++" 
    	```c++
		double real(Scalar const &s)
		```

!!! method "imag"

    Returns the imaginary part of the Scalar.
    === "C++" 
    	```c++
		double imag(Scalar const &s)
		```
		
!!! method "conj"

    Returns the complex conjugate of the Scalar.
    === "C++" 
    	```c++
		Scalar conj(Scalar const &s)
		```

!!! method "abs"

    Returns the absolute value of the Scalar.
    === "C++" 
    	```c++
		double abs(Scalar const &s)
		```

!!! method "isapprox"

    Compares two Scalars to whether they are approximately equal. 
	Return true if the following condition holds.

	$$ | a - b | < atol + rtol*|b| $$


    === "C++" 
    	```c++
		bool isapprox(Scalar const &a, Scalar const &b, double rtol = 1e-12, double atol = 1e-12)
		```

	**Parameters**
    
	| Name | Description        | Default |
    |------|--------------------|---------|
    | a    | first Scalar       |         |
    | b    | second Scalar      |         |
    | rtol | relative tolerance | 1e-12   |
	| atol | absolute tolerance | 0       |


!!! method "to_string"

	Converts the scalar into a string.
    === "C++" 
    	```c++
		std::string to_string(Scalar const &v);
		```
