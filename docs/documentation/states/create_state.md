---
title: Creating specific States
---

**Sources** [create_state.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/states/create_state.hpp), [create_state.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/states/create_state.cpp)

## product_state

Creates a filled product state. The local states are encoded using an integer. 

### Constructors

=== "Julia"
	```julia
	product_state(block::Block, local_states::Vector{Int64}; real::Bool=true)
	```
=== "C++"
	```c++
	State product_state(Block const &block, std::vector<int64_t> const &local_state, bool real = true);
	```
	
### Parameters

| Name         | Description                               |   |
|:-------------|:------------------------------------------|---|
| block        | block on which the state is defined       |   |
| local_states | local configurations of the product state |   |
| real         | flag whether real state is created        |   |

Each integer represents a certain local physical state which differs from block to block. Here is a summary for each one of the main 5 block types:

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

## random_state

Create a filled random state with normal $\mathcal{N}(0, 1)$ distributed coefficients.

=== "Julia"
	```julia
	random_state(block::Block; real::Bool=true, ncols::Int64=1, seed::Int64=42, normalized::Bool=true
	```
=== "C++"
	```c++
	State random_state(Block const &block, bool real = true, int64_t ncols = 1, int64_t seed = 42, bool normalized = true);
	```
	
### Parameters

| Name       | Description                                        |   |
|:-----------|:---------------------------------------------------|---|
| block      | block on which the state is defined                |   |
| real       | flag whether real state is created                 |   |
| ncols      | number of columns in the state                     |   |
| seed       | random seed determining the precise random numbers |   |
| normalized | flag whether the state is normalized               |   |

---

## zero_state

Create a filled state with all zero entries.

=== "Julia"
	```julia
	zero_state(block::Block; real::Bool=true, ncols::Int64=1)
	```
=== "C++"
	```c++
	State zero_state(Block const &block, bool real = true, int64_t ncols = 1);
	```
	
### Parameters

| Name  | Description                         |   |
|:------|:------------------------------------|---|
| block | block on which the state is defined |   |
| real  | flag whether real state is created  |   |
| ncols | number of columns in the state      |   |

---

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
	```julia
	--8<-- "examples/usage_examples/main.jl:create_state"
	```
=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:create_state"
	```



