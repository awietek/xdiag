---
title: OpSum
---

A sum of local operators acting on several lattice sites.

**Source** [opsum.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/operators/opsum.hpp)

## Constructor

=== "Julia"
	```julia
	OpSum()
	OpSum(ops::Vector{Op})
	```

=== "C++"	
	```c++
    OpSum() = default;
    OpSum(std::vector<Op> const &ops);
	```

| Parameter | Description                                                          |   |
|:----------|:---------------------------------------------------------------------|---|
| ops       | a vector of [Op](op.md) objects describing the operators summed over |   |

## Methods


!!! method "size"

	Returns the number of local [Op](op.md) operators.

	=== "Julia"
		```julia
		size(ops::OpSum)
		```

	=== "C++"	
		```c++
	    int64_t size() const;
		```

!!! method "defined"
	Returns bool whether a coupling (of type string) is defined

	=== "Julia"
		```julia
		defined(ops::OpSum, name::String)
		```

	=== "C++"	
		```c++
		bool defined(std::string name) const;
		```
		
!!! method "setindex! / operator[]"

	Sets a coupling given as a string to a certain numerical value or matrix

	=== "Julia"
		```julia
	    Base.setindex!(ops::OpSum, cpl, name::String)
		```

	=== "C++"	
		```c++
		Coupling &operator[](std::string name);
		```
		
!!! method "getindex / operator[]"

	Returns the value of a [Coupling](coupling.md) defined as a string. 

	=== "Julia"
		```julia
	    getindex(ops::OpSum, name::String)
		```

	=== "C++"	
		```c++
		Coupling const &operator[](std::string name) const;
		```

!!! method "couplings"

	Returns all the possible names of [Coupling](coupling.md) as a vector of strings

	=== "Julia"
		```julia
		couplings(ops::OpSum)
		```

	=== "C++"	
		```c++
		std::vector<std::string> couplings() const;
		```
		
!!! method "isreal"

	Returns whether or not the OpSum is real. This will throw an error if some [Coupling](coupling.md) are only defined as a string.

	=== "Julia"
		```julia
	    isreal(ops::OpSum)
		```

	=== "C++"	
		```c++
		bool isreal() const;
		```	
		
				
!!! method "isexplicit"

	Returns **false** if there exist a [Coupling](coupling.md) which is defined as a string, otherwise **true**.

	=== "Julia"
		```julia
		isexplicit(ops::OpSum)
		```

	=== "C++"	
		```c++
		bool isexplicit() const;
		```	
		
!!! method "operator +"

	Adds a single [Op](op.md) or a full OpSum

	=== "Julia"
		```julia
		+(ops::OpSum, op2::Op) = OpSum(ops.cxx_opsum + op2.cxx_op)
		+(ops::OpSum, ops2::OpSum) = OpSum(ops.cxx_opsum + ops2.cxx_opsum)
		```

	=== "C++"	
		```c++
        void operator+=(Op const &op);
		void operator+=(OpSum const &ops);
		OpSum operator+(Op const &op) const;
		OpSum operator+(OpSum const &ops) const;
		```	

## Usage Example

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:opsum"
	```

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:opsum"
	```

