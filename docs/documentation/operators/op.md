---
title: Op
---

A local operator acting on several lattice sites.

**Source** [op.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/operators/op.hpp)

## Constructors

=== "Julia"
	```julia
	Op(type::String, coupling::String, sites::Vector{Int64})
	Op(type::String, coupling::String, site::Int64)

	Op(type::String, coupling::Float64, sites::Vector{Int64})
	Op(type::String, coupling::Float64, site::Int64)

	Op(type::String, coupling::ComplexF64, sites::Vector{Int64})
	Op(type::String, coupling::ComplexF64, site::Int64)
	
	Op(type::String, coupling::Matrix{Float64}, sites::Vector{Int64})
	Op(type::String, coupling::Matrix{Float64}, site::Int64)

	Op(type::String, coupling::Matrix{ComplexF64}, sites::Vector{Int64})
	Op(type::String, coupling::Matrix{ComplexF64}, site::Int64)
	```

=== "C++"	
	```c++
    Op(std::string type, std::string coupling, std::vector<int64_t> const &sites)
    Op(std::string type, std::string coupling, int64_t site)
	
    Op(std::string type, double coupling, std::vector<int64_t> const &sites)
    Op(std::string type, double coupling, int64_t site)
	
    Op(std::string type, complex coupling, std::vector<int64_t> const &sites)
    Op(std::string type, complex coupling, int64_t site)
	
    Op(std::string type, arma::mat const &coupling, std::vector<int64_t> const &sites)
    Op(std::string type, arma::mat const &coupling, int64_t site)
	
    Op(std::string type, arma::cx_mat const &coupling, std::vector<int64_t> const &sites)
    Op(std::string type, arma::cx_mat const &coupling, int64_t site)
	```

| Parameter | Description                                                                |   |
|:----------|:---------------------------------------------------------------------------|---|
| type      | a string which denotes what kind of operator is represented                |   |
| coupling  | sets the coefficients neded to specify the coupling. Further details below |   |
| sites     | defines on which site(s) of the lattice the operator acts on.              |   |

The coupling can take on several types and allow some flexibility in defining operators.

| type                | Description                                                                         |   |
|:--------------------|:------------------------------------------------------------------------------------|---|
| string              | the coupling is represented as a string, e.g. "$t$" or "$J$" in a $t-J$ model |   |
| real/cplx number | the actual numerical value of the coupling                                          |   |
| real/cplx matrix | more generic interactions can be specified as matrices                              |   |


## Methods


!!! method "type"

	Returns the type of the operator

	=== "Julia"
		```julia
		type(op::Op)
		```

	=== "C++"	
		```c++
		std::string type() const;
		```

!!! method "coupling"
	Returns the coupling of the operator

	=== "Julia"
		```julia
		coupling(op::Op)
		```

	=== "C++"	
		```c++
		Coupling const &coupling() const;
		```
		
	This returns an object of type [Coupling](coupling.md), which can then be converted to an appropriate type.
		
!!! method "size"

	Returns how many sites the operator is defined on

	=== "Julia"
		```julia
		size(op::Op)
		```

	=== "C++"	
		```c++
		int64_t size() const;
		```
		
!!! method "getindex / operator[]"

	Returns the site with the given index.

	=== "Julia"
		```julia
	    getindex(op::Op, idx::Int64)
		```

	=== "C++"	
		```c++
		int64_t operator[](int64_t idx) const;
		```

!!! method "sites"

	Returns all the sites the operator is defined on

	=== "Julia"
		```julia
		sites(op::Op)
		```

	=== "C++"	
		```c++
		std::vector<int64_t> const &sites() const;
		```
		
!!! method "isreal"

	Returns whether or not the coupling is real. Throws an error if the coupling is given as a string, since then it cannot be determined whether the operator is real.

	=== "Julia"
		```julia
		isreal(op::Op)
		```

	=== "C++"	
		```c++
		bool isreal() const;
		```	
		
!!! method "ismatrix"

	Returns whether or not the coupling is defined as a matrix. Throws an error if the coupling is given as a string, since then it cannot be determined whether the operator is real.

	=== "Julia"
		```julia
		ismatrix(op::Op)
		```

	=== "C++"	
		```c++
		bool ismatrix() const;
		```	
		
				
!!! method "isexplicit"

	Returns **false** if the coupling is defined as a string, otherwise **true**

	=== "Julia"
		```julia
		isexplicit(op::Op)
		```

	=== "C++"	
		```c++
		bool isexplicit() const;
		```	

## Usage Example

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:op"
	```

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:op"
	```

