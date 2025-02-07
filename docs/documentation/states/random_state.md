---
title: RandomState
---

A random state with $\mathcal{N}(0, 1)$ normal distributed coefficients.

**Sources**<br> 
[random_state.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/states/random_state.hpp)<br>
[random_state.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/states/random_state.cpp)<br>
[random_state.jl](https://github.com/awietek/XDiag.jl/blob/main/src/states/random_state.jl)

---

## Constructors

=== "C++"	
	```c++
    RandomState(int64_t seed = 42, bool normalized = true);
	```
	
=== "Julia"
	```julia
	RandomState(seed::Int64 = 42, normalized::Bool = true)
	```

| Parameter  | Description                                          |   |
|:-----------|:-----------------------------------------------------|---|
| seed       | random seed determining which random numbers are put |   |
| normalized | flag whether the State is normalized                 |   |

---

## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:random_state"
	```

=== "Julia"
	```julia
	--8<-- "examples/usage_examples/main.jl:random_state"
	```

