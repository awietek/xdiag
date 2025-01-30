---
title: apply
---

Applies an operator given as an [Op](../operators/op.md) or [OpSum](../operators/opsum.md) to a [State](../states/state.md) $\vert w \rangle = \mathcal{O} \vert v\rangle$.

**Sources** [apply.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/states/apply.hpp), [apply.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/states/apply.hpp)

---

## Definition

An operator $\mathcal{O}$ can be applied to a state $\vert v\rangle$ in two ways:

1. Only the input state on which the operator acts is defined is given. The block of the output state is calculated and eventually created automatically.

	=== "C++"
		```c++
		State apply(Op const &op, State const &v);
		State apply(OpSum const &ops, State const &v);
		```
	=== "Julia"
		```julia
	    apply(op::Op, v::State)
	    apply(ops::OpSum, v::State)
		```

2. The output state is also handed as an argument which is overwritten. The compatibility of quantum numbers is checked. This way the output block is not created automatically and, thus, can be used to save computation time if the output block appears repeatedly in the computation.

	=== "C++"
		```c++
		void apply(Op const &op, State const &v, State &w);
		void apply(OpSum const &ops, State const &v, State &w);
		```
	=== "Julia"
		```julia
	    apply(op::Op, v::State, w::State)
	    apply(ops::OpSum, v::State, w::State)
		```

---

## Parameters

| Name     | Description                                                                      |   |
|:---------|:---------------------------------------------------------------------------------|---|
| ops / op | [OpSum](../operators/opsum.md) or [Op](../operators/op.md) defining the operator |   |
| v        | Input [State](../states/state.md) $\vert v\rangle  $                             |   |
| w        | Output [State](../states/state.md) $\vert w \rangle = O \vert v\rangle$          |   |

---

## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:apply"
	```

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:apply"
	```


