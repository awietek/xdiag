---
title: correlation_matrix
---

Computes the two-point correlations of a pair of single-site operators between all pairs of sites of a [State](state.md), returning the matrix

$$ C_{ij} = \langle \psi | \mathcal{O}^{(1)}_i \, \mathcal{O}^{(2)}_j | \psi \rangle. $$

For example, choosing the operator types `"Adag"` and `"A"` yields the single-particle density matrix.

**Sources:** [correlation_matrix.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/states/correlation_matrix.hpp) · [correlation_matrix.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/states/correlation_matrix.cpp)

---

## Definition

=== "C++"
	```c++
	arma::mat correlation_matrix(State const &state, std::string type1, std::string type2);
	arma::cx_mat correlation_matrixC(State const &state, std::string type1, std::string type2);
	```
=== "Julia"
	```julia
	correlation_matrix(state::State, type1::String, type2::String)
	```

The function `correlation_matrix` returns a real result and raises an error if the result is complex; in that case use `correlation_matrixC` in C++. In Julia, `correlation_matrix` returns a real or complex matrix as appropriate.

---

## Parameters

| Name  | Description                                                                                           |
|:------|:------------------------------------------------------------------------------------------------------|
| state | the [State](state.md) $|\psi\rangle$ the correlations are evaluated on (assumed normalized)            |
| type1 | the type of the first single-site [operator](../operators/operator_types.md) $\mathcal{O}^{(1)}$      |
| type2 | the type of the second single-site [operator](../operators/operator_types.md) $\mathcal{O}^{(2)}$     |

---

## Returns

An $N \times N$ matrix (with $N$ the number of sites), whose entry $(i, j)$ is $\langle \psi | \mathcal{O}^{(1)}_i \, \mathcal{O}^{(2)}_j | \psi \rangle$.

!!! info

	The state is assumed to be normalized; the raw matrix elements are returned without dividing by $\langle \psi | \psi \rangle$. The product of the two operators must keep the state within its block; a product that changes a conserved quantum number is reported as a block mismatch. On a symmetry-adapted block the operator is automatically symmetrized with the trivial irrep so that it acts within the block.

---

## Usage Example

=== "Julia"
	```julia
	block = Spinhalf(8)
	ops = OpSum()
	for i in 1:8
	    ops += "J" * Op("SdotS", [i, mod1(i + 1, 8)])
	end
	ops["J"] = 1.0
	e0, psi0 = eig0(ops, block)
	szsz = correlation_matrix(psi0, "Sz", "Sz")
	```
=== "C++"
	```c++
	auto block = Spinhalf(8);
	auto ops = OpSum();
	for (int i = 0; i < 8; ++i) {
	  ops += "J" * Op("SdotS", {i, (i + 1) % 8});
	}
	ops["J"] = 1.0;
	auto [e0, psi0] = eig0(ops, block);
	arma::mat szsz = correlation_matrix(psi0, "Sz", "Sz");
	```
