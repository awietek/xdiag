---
title: apply
---

Applies an operator to a state $\vert w \rangle = O \vert v\rangle$ or block of states $\left( \vert w_1 \rangle \dots \vert w_M \rangle \right) = O \left( \vert v_1 \rangle \dots \vert v_M \rangle \right)$. 

**Source** [apply.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/states/apply.hpp)

=== "Julia"
	```julia
	apply(op::Op, v::State, w::State, precision::Float64 = 1e-12)
	apply(ops::OpSum, v::State, w::State, precision::Float64 = 1e-12)
	```

=== "C++"
	```c++
    void apply(Op const &op, State const &v, State &w, double precision = 1e-12);
	void apply(OpSum const &ops, State const &v, State &w, double precision = 1e-12);
	```

The resulting state is handed as the third argument and is overwritten upon exit. 

## Parameters

| Name      | Description                                                                      |   |
|:----------|:---------------------------------------------------------------------------------|---|
| ops / op  | [OpSum](../operators/opsum.md) or [Op](../operators/op.md) defining the operator |   |
| v         | input state $\vert v\rangle  $                                                   |   |
| w         | output state $\vert w \rangle = O \vert v\rangle$                                |   |
| precision | precision with which checks for zero are performed (default $10^{-12}$)          |   |


## Usage Example


=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:apply"
	```

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:apply"
	```
	
