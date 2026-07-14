---
title: Measurements
---

# Measurements

Measurements in the form of expectation values of wavefunctions,

$$   \langle \mathcal{O}\rangle =  \langle \psi | \mathcal{O} | \psi \rangle,$$

can be evaluated using the [inner](../documentation/states/algebra.md#inner) function. For example, we compute a static spin correlation $\langle S_{0}^{z} S_{j}^{z}\rangle$ between site $0$ (resp. $1$ in Julia) and $j$.

=== "Julia"
	```julia
	--8<-- "examples/user_guide/main.jl:usage_guide_measu1"
	```
=== "C++"
	```c++
	--8<-- "examples/user_guide/main.cpp:usage_guide_measu1"
	```

Notice, that in Julia sites start counting from $1$, whereas in C++ sites are counted from $0$. Furthermore, if a complex wave function or operator is involved, the function [innerC](../documentation/states/algebra.md#inner) in C++ should be called, which returns a complex number. In Julia only [inner](../documentation/states/algebra.md#inner) is available whose return type is decided at runtime.

## Correlation functions

The [operator algebra](04-operators.md#the-operator-algebra) introduced earlier makes it straightforward to measure *arbitrary* correlation functions. Since operators can be multiplied, we can build any composite operator on the fly and hand it to [inner](../documentation/states/algebra.md#inner). For instance, the mixed spin correlation $\langle S^x_0 S^y_1 \rangle$ is obtained by multiplying the two single-site operators and measuring the resulting product.

=== "Julia"
	```julia
	--8<-- "examples/user_guide/main.jl:usage_guide_measu2"
	```
=== "C++"
	```c++
	--8<-- "examples/user_guide/main.cpp:usage_guide_measu2"
	```

This works for products on any number of sites and for any combination of operator types, allowing generic multi-point correlators to be assembled from elementary [Op](../documentation/operators/op.md) objects. Note that such correlators are frequently complex, so the [innerC](../documentation/states/algebra.md#inner) function is used in the C++ version.

## Local expectation values

A common task is to evaluate the expectation value of a single-site operator on *every* site of the lattice at once, for example the local magnetization $\langle S_i^z \rangle$ or the local density $\langle n_i \rangle$. This is conveniently done with the [expect](../documentation/states/expect.md) function, which takes the state and the type of a single-site operator and returns a vector whose $i$-th entry is $\langle \psi | \mathcal{O}_i | \psi \rangle$.

=== "Julia"
	```julia
	--8<-- "examples/user_guide/main.jl:usage_guide_expect"
	```
=== "C++"
	```c++
	--8<-- "examples/user_guide/main.cpp:usage_guide_expect"
	```

## Correlation matrices

Similarly, the two-point correlations of a pair of single-site operators between all pairs of sites can be computed in one call with the [correlation_matrix](../documentation/states/correlation_matrix.md) function. Given two operator types, it returns the matrix

$$ C_{ij} = \langle \psi | \mathcal{O}^{(1)}_i \, \mathcal{O}^{(2)}_j | \psi \rangle. $$

=== "Julia"
	```julia
	--8<-- "examples/user_guide/main.jl:usage_guide_corr"
	```
=== "C++"
	```c++
	--8<-- "examples/user_guide/main.cpp:usage_guide_corr"
	```

For instance, choosing the operator types `"Adag"` and `"A"` yields the single-particle density matrix. On a symmetry-adapted block, the operator is automatically symmetrized so that it acts within the block.

A few points are worth keeping in mind when using [expect](../documentation/states/expect.md) and [correlation_matrix](../documentation/states/correlation_matrix.md):

* The state is assumed to be **normalized**; the raw matrix elements are returned without dividing by $\langle \psi | \psi \rangle$.
* The operator (or the product of operators) must keep the state **within its block**. An operator that changes a conserved quantum number would make the expectation value vanish by symmetry and is instead reported as a block mismatch.
* Both functions return a real result and raise an error if the result is complex. In that case, the complex-valued variants [expectC](../documentation/states/expect.md) and [correlation_matrixC](../documentation/states/correlation_matrix.md) should be used in C++. In Julia, [expect](../documentation/states/expect.md) and [correlation_matrix](../documentation/states/correlation_matrix.md) return a real or complex result as appropriate.
