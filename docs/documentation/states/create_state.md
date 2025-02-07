---
title: Creating specific States
---

**Sources** [create_state.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/states/create_state.hpp), [create_state.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/states/create_state.cpp)

## product_state

Creates a filled product state.

=== "C++"
	```c++
	State product_state(Block const &block, std::vector<std::string> const &local_state, bool real = true);
	```
	
=== "Julia"
	```julia
	product_state(block::Block, local_states::Vector{String}, real::Bool=true)
	```


### Parameters

| Name         | Description                               |   |
|:-------------|:------------------------------------------|---|
| block        | block on which the state is defined       |   |
| local_states | local configurations of the product state |   |
| real         | flag whether real state is created        |   |

---

## random_state

Create a filled random state with normal $\mathcal{N}(0, 1)$ distributed coefficients.

=== "C++"
	```c++
	State random_state(Block const &block, bool real = true, int64_t seed = 42, bool normalized = true);
	```
	
=== "Julia"
	```julia
	random_state(block::Block, real::Bool=true, seed::Int64=42, normalized::Bool=true
	```


### Parameters

| Name       | Description                                        |   |
|:-----------|:---------------------------------------------------|---|
| block      | block on which the state is defined                |   |
| real       | flag whether real state is created                 |   |
| seed       | random seed determining the precise random numbers |   |
| normalized | flag whether the state is normalized               |   |

---

## zero_state

Create a filled state with all zero entries.

=== "C++"
	```c++
	State zero_state(Block const &block, bool real = true, int64_t n_cols = 1);
	```
	
=== "Julia"
	```julia
	zero_state(block::Block, real::Bool=true, n_col::Int64=1)
	```


### Parameters

| Name  | Description                         |   |
|:------|:------------------------------------|---|
| block | block on which the state is defined |   |
| real  | flag whether real state is created  |   |
| n_col | number of columns in the state      |   |

---

## zero

Set all coefficients of a given state to zero.

=== "C++"
	```c++
	void zero(State &state);
	```
	
=== "Julia"
	```julia
	zero(state::State)
	```
---

## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:create_state"
	```

=== "Julia"
	```julia
	--8<-- "examples/usage_examples/main.jl:create_state"
	```

