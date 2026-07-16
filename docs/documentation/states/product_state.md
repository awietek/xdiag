---
title: ProductState
---

A product state of local configurations on a given number of sites. The local states are encoded using an integer. Each integer represents a certain local physical state which differs from block to block. Here is a summary for each one of the main 5 block types:

#### Fermion
| Integer | Configuration     | Symbol |
|:--------|:------------------|:-------|
| `0`     | empty             | ○      |
| `1`     | occupied fermion  | ●      |

#### Boson
The local configuration of a site is simply the **occupation number**, i.e. an integer in $0, 1, \ldots, d-1$.

#### Spinhalf
| Integer | Configuration      | Symbol |
|:--------|:-------------------|:-------|
| `0`     | down-spin          | ↓      |
| `1`     | up-spin            | ↑      |

#### tJ
| Integer | Configuration          | Symbol |
|:--------|:-----------------------|:-------|
| `0`     | empty                  | ○      |
| `1`     | up-spin electron       | ↑      |
| `2`     | down-spin electron     | ↓      |

#### Electron
| Integer | Configuration          | Symbol |
|:--------|:-----------------------|:-------|
| `0`     | empty                  | ○      |
| `1`     | up-spin electron       | ↑      |
| `2`     | down-spin electron     | ↓      |
| `3`     | doubly occupied        | ⇅      |

For the fermionic blocks, [Fermion](../blocks/fermion.md),[tJ](../blocks/tJ.md), [Electron](../blocks/electron.md) a sign convention for normal ordering is chosen as explained in [normal ordering](../../user_guide/03-hilbert-spaces.md#normal-ordering-of-fermionic-blocks). 

**Sources:** [product_state.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/states/product_state.hpp) · [product_state.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/states/product_state.cpp)

## Constructors

=== "Julia"
	```julia
	ProductState(nsites::Int64)
	ProductState(local_states::Vector{Int64})
	```
=== "C++"	
	```c++
    ProductState(int64_t nsites);
	ProductState(std::vector<int64_t> const &local_states);
	```
	
| Parameter    | Description                                   |
|:-------------|:----------------------------------------------|
| nsites       | construct a product state on nsites           |
| local_states | the local configurations of the product state |

## Iteration

A ProductState can be iterated over, where at each iteration the string of the local configuration is retured. Here is an example:

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

#### nsites

Returns the number of sites of the product state

=== "Julia"
	```julia
	nsites(p::ProductState)
	```
=== "C++"	
	```c++
	int64_t nsites(ProductState const &p);
	```
	
#### size

Returns the number of sites of the product state. Same as "nsites".

=== "Julia"
	```julia
    size(state::ProductState)
	```
=== "C++"	
	```c++
	int64_t size(ProductState const &p);
	```
	
#### setindex! / operator[]

Sets the local configuration at the given site index to the given string.

=== "Julia"
	```julia
    setindex!(state::ProductState, local_state::String, idx::Int64)
	```
=== "C++"	
	```c++  
	std::string &operator[](int64_t i);
	```
	
#### getindex / operator[]

Returns the string of the local configuration at the given site index.
	
=== "Julia"
	```julia
    getindex(state::ProductState, idx::Int64)
	```
=== "C++"	
	```c++
	std::string const &operator[](int64_t i) const;
	```

#### push! / push_back

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
	```julia
	--8<-- "examples/usage_examples/main.jl:product_state"
	```
=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:product_state"
	```




