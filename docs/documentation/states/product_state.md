---
title: ProductState
---

A product state of local configurations

**Source** [product_state.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/states/product_state.hpp)

## Constructors

=== "Julia"
	```julia
	ProductState(nsites::Int64)
	ProductState(local_states::Vector{String})
	```

=== "C++"	
	```c++
    ProductState(int64_t nsites);
	ProductState(std::vector<std::string> const &local_states);
	```

| Parameter    | Description                                   |   |
|:-------------|:----------------------------------------------|---|
| nsites      | construct a product state on nsites          |   |
| local_states | the local configurations of the product state |   |

## Iteration

A Product state can be iterated over, where at each iteration the string of the local configuration is retured. Here is an example:

=== "Julia"
	```julia
	pstate = ProductState(["Up", "Dn", "Emp", "UpDn"])
	for s in pstate
		@show s
	end
	```

=== "C++"	
	```c++
    auto pstate = ProductState({"Up", "Dn", "Emp", "UpDn"});
	for (auto s : pstate) {
		Log("{}", s);
	}
	```
	
## Methods

!!! method "nsites"

	Returns the number of sites of the product state

	=== "Julia"
		```julia
		nsites(state::ProductState)
		```

	=== "C++"	
		```c++
		int64_t nsites() const
		```

!!! method "size"

	Returns the number of sites of the product state. Same as "nsites".

	=== "Julia"
		```julia
	    size(state::ProductState)
		```

	=== "C++"	
		```c++
		int64_t size() const
		```

!!! method "setindex! / operator[]"

	Sets the local configuration at the given site index to the given string.

	=== "Julia"
		```julia
	    setindex!(state::ProductState, local_state::String, idx::Int64)
		```

	=== "C++"	
		```c++  
		std::string &operator[](int64_t i);
		```
		
!!! method "getindex / operator[]"

	Returns the string of the local configuration at the given site index.

	=== "Julia"
		```julia
	    getindex(state::ProductState, idx::Int64)
		```

	=== "C++"	
		```c++
		std::string const &operator[](int64_t i) const;
		```

!!! method "push! / push_back"

	Adds a local configuration add the end of the product state.

	=== "Julia"
		```julia
	    push!(state::ProductState, local_state::String
		```

	=== "C++"	
		```c++
		void push_back(std::string l);
		```

## Usage Example

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:product_state"
	```

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:product_state"
	```

