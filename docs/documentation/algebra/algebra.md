---
title: Algebra routines 
---

Several basic algebraic operations for states and operators.

**Sources**<br>
[algebra.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/algebra/algebra.hpp)<br>
[algebra.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/algebra/algebra.cpp)<br>
[algebra.jl](https://github.com/awietek/XDiag.jl/blob/main/src/algebra/algebra.jl)

---

## dot

Computes the dot product $\langle v \vert w \rangle$ between two states $\vert v \rangle$ and $\vert w \rangle$. In C++, please use the dotC function if one of the two states is expected to be complex.

=== "C++"
	```c++
    double dot(State const &v, State const &w);
	complex dotC(State const &v, State const &w);
	```
	
=== "Julia"
	```julia
	dot(v::State, w::State)
	```

---

## inner

Computes the expectation value $\langle v | O |v \rangle$ of an operator $O$ and a state $|v\rangle$. The operator can either be an [Op]("../operators/op.md") or an [OpSum]("../operators/opsum.md") object. In C++, please use the innerC function if either the operator or the state are complex.

=== "C++"
	```c++
	double inner(Op const &op, State const &v);
	double inner(OpSum const &ops, State const &v);
	complex innerC(Op const &op, State const &v);
	complex innerC(OpSum const &ops, State const &v);
	```

=== "Julia"
	```julia
	inner(op::Op, v::State)
	inner(ops::OpSum, v::State)
	```

---

## norm

Computes the 2-norm $\parallel |v \rangle \parallel_2$ of a state $|v \rangle$ defined as

$$ \parallel |v \rangle \parallel_2 = \sum_n |\langle n | v \rangle |^2, $$

where $\{ |n\rangle \}$ denotes an orthonormal basis of the block.

=== "C++"
	```c++
    double norm(State const &v);
	```

=== "Julia"
	```julia
	norm(state::State)::Float64
	```
	
---

## norm1

Computes the 1-norm $\parallel |v \rangle \parallel_1$ of a state $|v \rangle$ defined as

$$ \parallel |v \rangle \parallel_1 = \sum_n |\langle n | v \rangle |, $$

where $\{ |n\rangle \}$ denotes an orthonormal basis of the block.


=== "C++"
	```c++
    double norm1(State const &v);
	```
	
=== "Julia"
	```julia
	norm1(state::State)::Float64
	```
	
---

## norminf


Computes the $\infty$-norm $\parallel |v \rangle \parallel_\infty$ of a state $|v \rangle$ defined as

$$ \parallel |v \rangle \parallel_\infty = \max_n |\langle n | v \rangle |, $$

where $\{ |n\rangle \}$ denotes an orthonormal basis of the block.

=== "C++"
	```c++
    double norminf(State const &v);
	```
	
=== "Julia"
	```julia
	norminf(state::State)::Float64
	```
	
---

## Usage Examples

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:algebra"
	```

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:algebra"
	```
	
