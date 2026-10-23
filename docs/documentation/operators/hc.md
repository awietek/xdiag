---
title: hc, adjoint
---

**Sources:** [hc.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/operators/hc.hpp) · [hc.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/operators/hc.cpp)

Returns the hermitian conjugate $\mathcal{O}^\dagger$ of an operator $\mathcal{O}$ represented by an [Op](op.md) or [OpSum](opsum.md) object, such that on every block
$$ \mathrm{matrix}(\mathrm{hc}(\mathcal{O})) = \mathrm{matrix}(\mathcal{O})^\dagger. $$
`adjoint` is the same operation under the name used in the language of $*$-algebras; in Julia it also provides the postfix `'` operator, so `ops'` is equivalent to `hc(ops)`. Please note the details when conjugating complex couplings, outlined in [OpSum # Complex couplings](opsum.md#complex-couplings).

## Definition

=== "Julia"
	```julia
	hc(op::Op)::OpSum
	hc(ops::OpSum)::OpSum
	Base.adjoint(op::Op)::OpSum
	Base.adjoint(ops::OpSum)::OpSum
	```
=== "C++"
	```c++
	OpSum hc(Op const &op)
	OpSum hc(OpSum const &ops)
	OpSum adjoint(Op const &op)
	OpSum adjoint(OpSum const &ops)
	```

Hermitian conjugation complex-conjugates the couplings and reverses the order of products,
$$ \mathrm{hc}(c\,\mathcal{A}\mathcal{B}) = c^*\, \mathcal{B}^\dagger \mathcal{A}^\dagger. $$

The hermitian conjugate of an [Op](op.md) need not necessarily be of type [Op](op.md) as well, since for certain operator types, a numerical prefactor can be introduced (e.g. `HopAsym` or `ExchangeAsym` introduce a `-1` upon hermitian conjugation, see [operator types](operator_types.md)) which then generically requires representing the object as an [OpSum](opsum.md). The elementary [operator types](operator_types.md) are mapped onto their hermitian-conjugate partners, e.g. $\mathrm{hc}(S^+) = S^-$, $\mathrm{hc}(c_\uparrow) = c^\dagger_\uparrow$, while hermitian types such as $S^z$, $S^y$ or `SdotS` are unchanged. For an explicit [`Matrix`](operator_types.md) type, the conjugate transpose of the given matrix is taken.

## Usage Example

=== "Julia"
	```julia
	--8<-- "examples/usage_examples/main.jl:hc"
	```
=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:hc"
	```
