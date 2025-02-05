---
title: hc
---

Returns the hermitian conjugate $\mathcal{O}^\dagger$ of an operator $\mathcal{O}$ represented by an [Op](op.md) or [OpSum](opsum.md) object. Please note the details when conjugating complex couplings, outlined in [OpSum # Complex couplings](opsum.md#complex-couplings).

**Sources**<br>
[hc.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/operators/logic/hc.hpp)<br>
[hc.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/operators/logic/hc.cpp)<br>
[hc.jl](https://github.com/awietek/XDiag.jl/blob/main/xdiag/operators/logic/hc.jl)

---

## Definition

=== "C++"
	```c++
	Op hc(Op const &op)
	OpSum hc(OpSum const &ops)
	```
=== "Julia"
	```julia
	hc(op::OpSum)
	hc(ops::OpSum)
	```
	
---

## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:hc"
	```

=== "Julia"
	```julia
	--8<-- "examples/usage_examples/main.jl:hc"
	```
