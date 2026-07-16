---
title: Operators
---

# Operators

Besides Hilbert spaces, the second key objects in quantum mechanics are operators. In a many-body setting, we consider operators of the form,
$$
	O = \sum_{A\subseteq \mathcal{L}} c_A O_A,
$$
where $O_A$ denotes a local operator acting on sites $A=\{a_1, \ldots, a_{l_A}\}$ and $\mathcal{L}$ denotes the lattice and $c_{A}$ are coupling constants. In the case of the Heisenberg model, we would thus have $\mathcal{O}_{A} = \mathbf{S}_i\cdot\mathbf{S}_j$ and $c_A = J$. In XDiag, the local operators are represented via an [Op](../documentation/operators/op.md) object while the sum is represented by an [OpSum](../documentation/operators/opsum.md) object. These values of the coupling constants $c_A$ can either be a real or complex number, or a string that later needs to be replaced. The Hamiltonian of a spin $S=1 / 2$ Heisenberg chain is created in the following way:

=== "Julia"
	```julia
	--8<-- "examples/user_guide/main.jl:usage_guide_op1"
	```
=== "C++"
	```c++
	--8<-- "examples/user_guide/main.cpp:usage_guide_op1"
	```

We first create an empty [OpSum](../documentation/operators/opsum.md) and then add additional terms to it. The first part of the product denotes the coupling constant, here given as a string. Alternatively, one could have directly used real / complex numbers here. The second part of the product is a single [Op](../documentation/operators/op.md) object. It is created with two inputs:

* The type, here `SdotS` denoting an operator of the form $\mathbf{S}_{i} \cdot \mathbf{S}_{j}$. XDiag features a wide variety of different operator types.
* An array defining which site the operator lives on. Notice, that in Julia we start counting the sites from 1, while in C++ we start counting the sites from 0.

Which operator types are available, and on which blocks they may be used, is documented in two places:

* On the documentation page of each **block type** (for example [Spinhalf](../documentation/blocks/spinhalf.md), [tJ](../documentation/blocks/tJ.md), [Electron](../documentation/blocks/electron.md), [Boson](../documentation/blocks/boson.md), or [Fermion](../documentation/blocks/fermion.md)), a table lists all operator types that can be applied to that particular block.
* A central register of every operator type, together with its definition, its required number of sites, and the blocks it is defined for, is collected on the [operator types](../documentation/operators/operator_types.md) page.

## Coupling constants

The coefficient $c_A$ multiplying an operator can be given in two ways. It can be a plain **real or complex number**, which is the most direct choice. Alternatively, it can be a **named coupling** given as a string, such as `"J"` in the example above. Named couplings are convenient because they let us define the structure of a model once and only fix the numerical values later. This is especially useful when the same Hamiltonian is diagonalized repeatedly for different parameter values, or when parameters are read from an input file.

The value of a named coupling is assigned using the `[]` operator, and the currently assigned values can be queried the same way. Once all named couplings have been assigned a value, the `plain()` function returns a copy of the [OpSum](../documentation/operators/opsum.md) in which every named coupling has been substituted by its numerical value.

=== "Julia"
	```julia
	--8<-- "examples/user_guide/main.jl:usage_guide_op2"
	```
=== "C++"
	```c++
	--8<-- "examples/user_guide/main.cpp:usage_guide_op2"
	```

Most computational routines call `plain()` internally, so in practice it is rarely necessary to invoke it by hand. An [OpSum](../documentation/operators/opsum.md) that still contains unresolved named couplings cannot be turned into a matrix or applied to a state.

## The operator algebra

Beyond forming linear combinations, XDiag now supports the full set of algebraic operations on operators. In addition to **addition**, **subtraction**, and **multiplication by a scalar**, two [OpSum](../documentation/operators/opsum.md) objects can be **multiplied** with one another. This product is the ordinary composition of operators and is distributed over the sum,

$$ \Big(\sum_i c_i M_i\Big)\Big(\sum_j d_j N_j\Big) = \sum_{i,j} (c_i d_j)\, M_i N_j. $$

Together with the vector space operations, this multiplication turns the set of [OpSum](../documentation/operators/opsum.md) objects into an **algebra**. Note that the product is generally **non-commutative**, reflecting the fact that quantum mechanical operators need not commute. A simple example is the product of two single-site spin operators,

=== "Julia"
	```julia
	--8<-- "examples/user_guide/main.jl:usage_guide_op3"
	```
=== "C++"
	```c++
	--8<-- "examples/user_guide/main.cpp:usage_guide_op3"
	```

This makes it possible to build composite operators, such as correlation functions or powers of a Hamiltonian, directly from elementary building blocks.

## Hermitian conjugation

The algebra is complemented by an operation of **hermitian conjugation**, computed with the [hc](../documentation/operators/hc.md) function. It returns the hermitian conjugate $\mathcal{O}^\dagger$ of an [Op](../documentation/operators/op.md) or [OpSum](../documentation/operators/opsum.md), complex-conjugating the coefficients and conjugating each elementary operator (for example $\mathrm{hc}(S^+) = S^-$).

=== "Julia"
	```julia
	--8<-- "examples/user_guide/main.jl:usage_guide_op4"
	```
=== "C++"
	```c++
	--8<-- "examples/user_guide/main.cpp:usage_guide_op4"
	```

Hermitian conjugation is an *involution*: applying it twice returns the original operator, $(\mathcal{O}^\dagger)^\dagger = \mathcal{O}$, and it reverses the order of products, $(\mathcal{A}\mathcal{B})^\dagger = \mathcal{B}^\dagger \mathcal{A}^\dagger$. An algebra equipped with such an operation is called an **involutive algebra** (or $*$-algebra). The [OpSum](../documentation/operators/opsum.md) objects of XDiag therefore realize precisely the mathematical structure that operators in quantum mechanics are expected to have, and `hc` is the natural tool for building hermitian Hamiltonians, for instance by symmetrizing a non-hermitian term as `A + hc(A)`.
