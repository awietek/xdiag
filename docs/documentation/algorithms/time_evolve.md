---
title: time_evolve
---

Computes the real-time evolution, 

$$\vert \psi(t) \rangle = e^{-iHt} \vert \psi_0\rangle,$$ 

of a [State](../states/state.md) $\vert \psi_0 \rangle$ and a Hermitian operator $H$ using an iterative algorithm. 

**Sources** [time_evolve.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/algorithms/time_evolve.hpp), [time_evolve.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/algorithms/time_evolve.cpp)

## Definition

The method is provided in two variants:

1. Returning a new state while the input state remains untouched. This variant is safe to use and simple to code.

	=== "C++"

		```c++
		State time_evolve(OpSum const &H, State psi0, double time,
                          double precision = 1e-12,
                          std::string algorithm = "lanczos");
		```

2. An *inplace* variant `time_evolve_inplace`, where the input state is overwritten and contains the time evolved state upon exit. This version is more memory efficient than `time_evolve`.

	=== "C++"

		```c++
		void time_evolve_inplace(OpSum const &H, State &psi0, double time,
                                 double precision = 1e-12,
                                 std::string algorithm = "lanczos");
		```

## Parameters

| Name      | Description                                                                           | Default   |
|:----------|:--------------------------------------------------------------------------------------|-----------|
| H         | [OpSum](../operators/opsum.md) defining the hermitian operator $H$ for time evolution |           |
| psi0      | initial [State](../states/state.md) $\vert \psi_0 \rangle$ of the time evolution      |           |
| time      | time $t$ until which the state is evolved                                             |           |
| precision | accuracy of the computed time evolved state $\vert \psi(t) \rangle$                    | 1e-12     |
| algorithm | iterative algorithm which is used, one of `lanczos` or `expokit`                      | `lanczos` |

The `algorithm` parameter decised which backend is run. If `lanczos` is chosen, the [evolve_lanczos](evolve_lanczos.md) routine is called with the standard arguments. Alternatively, `expokit` chooses the [time_evolve_expokit](time_evolve_expokit.md) routine. For a detailed documentation of the algorithms we refer to the [evolve_lanczos](evolve_lanczos.md) and [time_evolve_expokit](time_evolve_expokit.md) pages. Broadly speaking, the `expokit` can yield higher precision states at arbitrarily long times at the cost of increased memory and computing time. In practice, we recommend analysing the effect of the `precision` parameters on the time evolution series obtained in both cases. 

## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:time_evolve"
	```
	
