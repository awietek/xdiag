---
title: fill
---

Fills a [State](state.md) with a given model state, e.g. a [ProductState](product_state.md) or a [RandomState](random_state.md).

**Source** [fill.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/states/fill.hpp)

## Definition

=== "Julia"
	```julia
	fill(state::State, pstate::ProductState, ncol::Int64 = 1)
	fill(state::State, rstate::RandomState, ncol::Int64 = 1)
	```

=== "C++"
	```c++
	void fill(State &state, ProductState const &pstate, int64_t ncol = 0);
	void fill(State &state, RandomState const &rstate, int64_t ncol = 0);
	```
## Parameters

| Name   | Description                                                                     |   |
|:-------|:--------------------------------------------------------------------------------|---|
| state  | [State](state.md) object to be filled                                           |   |
| pstate | [ProductState](product_state.md) object                                         |   |
| rstate | [RandomState](random_state.md) object                                           |   |
| ncol   | integer deciding which column of the State is filled (default: 1/0 (Julia/C++)) |   |



## Usage Example

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:fill"
	```

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:fill"
	```
	
