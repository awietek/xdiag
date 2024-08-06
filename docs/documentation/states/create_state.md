---
title: Functions to create and modify states
---

**Source** [create_state.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/states/create_state.hpp)

## product

Creates a filled product state.

=== "Julia"
	```julia
	product(block::Block, local_states::Vector{String}, real::Bool=true)
	```

=== "C++"
	```c++
	State product(Block const &block, std::vector<std::string> const &local_state, bool real = true);
	```
### Parameters

| Name         | Description                               |   |
|:-------------|:------------------------------------------|---|
| block        | block on which the state is defined       |   |
| local_states | local configurations of the product state |   |
| real         | flag whether real state is created        |   |

## rand

Create a filled random state with normal distributed coefficients.

=== "Julia"
	```julia
	rand(block::Block, real::Bool=true, seed::Int64=42, normalized::Bool=true
	```

=== "C++"
	```c++
	State rand(Block const &block, bool real = true, int64_t seed = 42, bool normalized = true);
	```
### Parameters

| Name       | Description                                        |   |
|:-----------|:---------------------------------------------------|---|
| block      | block on which the state is defined                |   |
| real       | flag whether real state is created                 |   |
| seed       | random seed determining the precise random numbers |   |
| normalized | flag whether the state is normalized               |   |

## zeros

Create a filled state with all zero entries.

=== "Julia"
	```julia
	zeros(block::Block, real::Bool=true, n_col::Int64=1)
	```

=== "C++"
	```c++
	State zeros(Block const &block, bool real = true, int64_t n_cols = 1);
	```
### Parameters

| Name  | Description                         |   |
|:------|:------------------------------------|---|
| block | block on which the state is defined |   |
| real  | flag whether real state is created  |   |
| n_col | number of columns in the state      |   |


## zero

Set all coefficients of a given state to zero.

=== "Julia"
	```julia
	zero(state::State)
	```

=== "C++"
	```c++
	void zero(State &state);
	```


## Usage Example

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:create_state"
	```

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:create_state"
	```
	
