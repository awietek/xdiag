---
title: Basic algebra routines 
---

**Source** [algebra.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/states/algebra.hpp)

## norm

Computes the 2-norm of a state

=== "Julia"
	```julia
	norm(state::State)
	```

=== "C++"
	```c++
    double norm(State const &v);
	```

## norm1

Computes the 1-norm of a state

=== "Julia"
	```julia
	norm1(state::State)
	```

=== "C++"
	```c++
    double norm1(State const &v);
	```

## norminf

Computes the $\infty$-norm of a state

=== "Julia"
	```julia
	norminf(state::State)
	```

=== "C++"
	```c++
    double norminf(State const &v);
	```


## dot

Computes the dot product between two states. In C++, please use the dotC function if one of the two states is expected to be complex.

=== "Julia"
	```julia
	dot(v::State, w::State)
	```

=== "C++"
	```c++
    double dot(State const &v, State const &w);
	complex dotC(State const &v, State const &w);
	```

## inner

Computes the expectation value $\langle v | O |v \rangle$ of an operator $O$ and a state $|v\rangle$. The operator can either be an [Op]("../operators/op.md") or an [OpSum]("../operators/opsum.md") object. In C++, please use the innerC function if either the operator or the state are complex.

=== "Julia"
	```julia
	inner(ops::OpSum, v::State)
	inner(op::Op, v::State)
	```

=== "C++"
	```c++
	double inner(OpSum const &ops, State const &v);
	double inner(Op const &op, State const &v);
	complex innerC(OpSum const &ops, State const &v);
	complex innerC(Op const &op, State const &v);
	```

## Usage Examples

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:algebra"
	```

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:algebra"
	```
	
