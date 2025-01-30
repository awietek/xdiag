---
title: ProductState
---

A product state of local configurations

**Sources** [product_state.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/states/product_state.hpp), [product_state.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/states/product_state.cpp)

--- 

## Constructors

=== "C++"	
	```c++
    ProductState(int64_t nsites);
	ProductState(std::vector<std::string> const &local_states);
	```
	
=== "Julia"
	```julia
	ProductState(nsites::Int64)
	ProductState(local_states::Vector{String})
	```



| Parameter    | Description                                   |
|:-------------|:----------------------------------------------|
| nsites       | construct a product state on nsites           |
| local_states | the local configurations of the product state |


---

## Iteration

A ProductState can be iterated over, where at each iteration the string of the local configuration is retured. Here is an example:

=== "C++"	
	```c++
    auto pstate = ProductState({"Up", "Dn", "Emp", "UpDn"});
	for (auto s : pstate) {
		Log("{}", s);
	}
	```
	
=== "Julia"
	```julia
	pstate = ProductState(["Up", "Dn", "Emp", "UpDn"])
	for s in pstate
		@show s
	end
	```

---
	
## Methods

#### nsites

Returns the number of sites of the product state

=== "C++"	
	```c++
	int64_t nsites(ProductState const &p);
	```
	
=== "Julia"
	```julia
	nsites(p::ProductState)
	```

---

#### size

Returns the number of sites of the product state. Same as "nsites".

=== "C++"	
	```c++
	int64_t size(ProductState const &p);
	```
	
=== "Julia"
	```julia
    size(state::ProductState)
	```

---

#### setindex! / operator[]

Sets the local configuration at the given site index to the given string.

=== "C++"	
	```c++  
	std::string &operator[](int64_t i);
	```
	
=== "Julia"
	```julia
    setindex!(state::ProductState, local_state::String, idx::Int64)
	```

---
	
#### getindex / operator[]

Returns the string of the local configuration at the given site index.

=== "C++"	
	```c++
	std::string const &operator[](int64_t i) const;
	```
	
=== "Julia"
	```julia
    getindex(state::ProductState, idx::Int64)
	```
---

#### push! / push_back

Adds a local configuration add the end of the product state.

=== "C++"	
	```c++
	void push_back(std::string l);
	```

=== "Julia"
	```julia
    push!(state::ProductState, local_state::String
	```

---

## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:product_state"
	```

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:product_state"
	```


