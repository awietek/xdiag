---
title: apply
---

Applies an operator given as an [Op](../operators/op.md) or [OpSum](../operators/opsum.md) to a [State](../states/state.md) $\vert w \rangle = \mathcal{O} \vert v\rangle$.

**Sources**<br>
[apply.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/algebra/apply.hpp)<br> 
[apply.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/algebra/apply.hpp)<br>
[apply.jl](https://github.com/awietek/XDiag.jl/blob/main/src/algebra/apply.jl)

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

3. A low-level routine where given an [OpSum](../operators/opsum.md) with an input and output block, the operator is applied to a raw vector. This version is useful when using XDiag in conjuction with third-party libraries. 

	=== "C++"
		```c++
		void apply(OpSum const &ops, 
		           Block const &block_in, arma::vec const &v, 
				   Block const &block_out, arma::vec &w);
		void apply(OpSum const &ops, 
			       Block const &block_in, arma::cx_vec const &v, 
				   Block const &block_out, arma::cx_vec &w);
		void apply(OpSum const &ops, 
		           Block const &block_in, arma::mat const &V, 
				   Block const &block_out, arma::mat &W);
		void apply(OpSum const &ops, 
		           Block const &block_in, arma::cx_mat const &V, 
				   Block const &block_out, arma::cx_mat &W);
		```
	=== "Julia"
		```julia
	    apply(ops::OpSum,
              block_in::Block, v::Vector{Float64},
              block_out::Block, w::Vector{Float64})
	    apply(ops::OpSum,
              block_in::Block, v::Vector{ComplexF64},
              block_out::Block, w::Vector{ComplexF64})
		apply(ops::OpSum,
               block_in::Block, V::Matrix{Float64},
               block_out::Block, W::Matrix{Float64})
		apply(ops::OpSum,
               block_in::Block, V::Matrix{ComplexF64},
               block_out::Block, W::Matrix{ComplexF64})
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


