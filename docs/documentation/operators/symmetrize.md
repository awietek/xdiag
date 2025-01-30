---
title: symmetrize
---
Symmetrizes an operator with respect to a [PermutationGroup](../symmetries/permutation_group.md) or a [Representation](../symmetries/representation.md).

**Sources** [symmetrize.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/operators/symmetrize.hpp), [symmetrize.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/operators/symmetrize.cpp)

---

Symmetrization in this context means the following. In general, we are given an [OpSum](../operators/opsum.md) of the form,

$$ O = \sum_{A\subseteq \mathcal{L}} O_A,$$

where $O_A$ denotes a local operator acting on sites $A=\{a_1, \ldots, a_{l_A}\}$ and $L$ denotes the lattice.
A [PermutationGroup](../symmetries/permutation_group.md) $\mathcal{G}$ is defined through its
permutations $\pi_1, \ldots, \pi_M$. The symmetrized operator returned by this function is then 

$$ O^\mathcal{G} = \frac{1}{M}\sum_{A\subseteq \mathcal{L}} \sum_{\pi \in \mathcal{G}}  O_{\pi(A)},$$

where $\pi(A) = \{\pi(a_1), \ldots,\pi(a_{l_A})\}$ denotes the permutated set of sites of the local operator $O_A$. If a [Representation](../symmetries/representation.md) called $\rho$ is given in addition, the following operator is constructed,

$$ O^\mathcal{G, \rho} = \frac{1}{M}\sum_{A\subseteq \mathcal{L}} \sum_{\pi \in \mathcal{G}} \chi_\rho(\pi) O_{\pi(A)},$$

where $\chi_\rho(\pi)$ denotes the characters of the representation $\rho$. This routine is useful to evaluate observables in symmetrized blocks.

---

## Definition

=== "C++"
	```c++
    OpSum symmetrize(Op const &op, PermutationGroup const &group);
	OpSum symmetrize(Op const &op, Representation const &irrep);
	OpSum symmetrize(OpSum const &ops, PermutationGroup const &group);
	OpSum symmetrize(OpSum const &ops, Representation const &irrep);
	```
	
=== "Julia"
	```julia
	symmetrize(op::Op, group::PermutationGroup)
	symmetrize(op::Op, irrep::Representation)
	symmetrize(ops::OpSum, group::PermutationGroup)
	symmetrize(ops::OpSum, irrep::Representation)
	```
	
---

## Parameters

| Name     | Description                                                                                        |   |
|:---------|:---------------------------------------------------------------------------------------------------|---|
| ops / op | [OpSum](../operators/opsum.md) or [Op](../operators/op.md) defining the operator to be symmetrized |   |
| group    | [PermutationGroup](../symmetries/permutation_group.md) defining the permutation symmetries         |   |
| irrep    | Irreducible [Representation](../symmetries/representation.md)  of the symmetry group               |   |

---

## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:symmetrize"
	```

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:symmetrize"
	```


