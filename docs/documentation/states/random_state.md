---
title: RandomState
---

A random state with normal distributed coefficients

**Source** [random_state.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/states/random_state.hpp)

## Constructors

=== "Julia"
	```julia
	RandomState(seed::Int64 = 42, normalized::Bool = true)
	```

=== "C++"	
	```c++
    RandomState(int64_t seed = 42, bool normalized = true);
	```

| Parameter  | Description                                          |   |
|:-----------|:-----------------------------------------------------|---|
| seed       | random seed determining which random numbers are put |   |
| normalized | flag whether the State is normalized                 |   |


	
## Methods

!!! method "seed"

	Returns the seed of the random state

	=== "Julia"
		```julia
		seed(state::RandomState)
		```

	=== "C++"	
		```c++
		int64_t seed() const;
		```

!!! method "size"

	Returns whether the state is normalized.

	=== "Julia"
		```julia
	    normalized(state::RandomState)
		```

	=== "C++"	
		```c++
		bool normalized() const;
		```


## Usage Example

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:random_state"
	```

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:random_state"
	```

