---
title: expect
---

Computes the expectation value of a single-site operator on *every* site of a [State](state.md), returning a vector whose $i$-th entry is $\langle \psi | \mathcal{O}_i | \psi \rangle$. This is convenient for evaluating quantities such as the local magnetization $\langle S^z_i \rangle$ or the local density $\langle n_i \rangle$ in a single call.

**Sources:** [expect.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/states/expect.hpp) · [expect.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/states/expect.cpp)

## Definition

=== "Julia"
	```julia
	expect(state::State, type::String)
	```
=== "C++"
	```c++
	arma::vec expect(State const &state, std::string type);
	arma::cx_vec expectC(State const &state, std::string type);
	```

The function `expect` returns a real result and raises an error if the result is complex; in that case use `expectC` in C++. In Julia, `expect` returns a real or complex vector as appropriate.

## Parameters

| Name  | Description                                                                                         |   |
|:------|:----------------------------------------------------------------------------------------------------|---|
| state | the [State](state.md) $\lvert\psi\rangle$ the expectation value is evaluated on (assumed normalized)    |   |
| type  | the type of a single-site [operator](../operators/operator_types.md) (e.g. `"Sz"`, `"N"`, `"Ntot"`) |   |

## Returns

A vector of length `nsites`, whose $i$-th entry is $\langle \psi | \mathcal{O}_i | \psi \rangle$ with $\mathcal{O}_i$ the operator of the given `type` acting on site $i$.

!!! info

	The state is assumed to be normalized; the raw matrix elements are returned without dividing by $\langle \psi | \psi \rangle$. The operator must keep the state within its block; an operator that changes a conserved quantum number is reported as a block mismatch.

## Usage Example

=== "Julia"
	```julia
	--8<-- "examples/usage_examples/main.jl:expect"
	```
=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:expect"
	```
