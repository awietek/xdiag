---
title: hc
---

Returns the hermitian conjugate $\mathcal{O}^\dagger$ of an operator $\mathcal{O}$ represented by an [Op](op.md) or [OpSum](opsum.md) object. Please note the details when conjugating complex couplings, outlined in [OpSum # Complex couplings](opsum.md#complex-couplings).

**Sources:** [hc.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/operators/hc.hpp) · [hc.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/operators/hc.cpp)

## Definition

=== "Julia"
	```julia
	hc(op::Op)::OpSum
	hc(ops::OpSum)::OpSum
	```
=== "C++"
	```c++
	OpSum hc(Op const &op)
	OpSum hc(OpSum const &ops)
	```

The hermitian conjugate of an [Op](op.md) need not necessarily be of type [Op](op.md) as well, since for certain operator types, a numerical prefactor can be introduced (e.g. `HopAsym` or `ExchangeAsym` introduce a `-1` upon hermitian conjugation, see [operator types](operator_types.md)) which then generically requires representing the object as an [OpSum](opsum.md).


## Usage Example

=== "Julia"
	```julia
	--8<-- "examples/usage_examples/main.jl:hc"
	```
=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:hc"
	```
