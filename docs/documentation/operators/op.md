---
title: Op
---

Object describing a single linear operator acting on a Hilbert space.

**Sources**<br>
[op.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/operators/op.hpp)<br>
[op.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/operators/op.cpp)<br>
[op.jl](https://github.com/awietek/XDiag.jl/blob/main/src/operators/op.jl)

---

Every operator is defined by up to three paramaters:

1. The **type** of the operator. This is a string argument which determines what kind of operator is represented. A detailed overview of the available types can be found at [Operator types](operator_types.md)

2. The **sites** of the operator. This defines which physical sites (or orbitals) the operator acts upon. While most operator types require sites, there are also operator types (e.g. `HubbardU`) which do not need to define sites.

3. For special interactions, it can be necessary to additionally specify a numerical **matrix**, which can either be real or complex. An example is the operator type `Matrix` defining generic spin interactions.

---

## Constructors

=== "C++"	
	```c++
	Op(std::string type);
	Op(std::string type, int64_t site);
	Op(std::string type, std::vector<int64_t> const &sites);
	Op(std::string type, int64_t site, arma::mat const &matrix);
	Op(std::string type, std::vector<int64_t> const &sites, arma::mat const &matrix);
	Op(std::string type, int64_t site, arma::cx_mat const &matrix);
	Op(std::string type, std::vector<int64_t> const &sites, arma::cx_mat const &matrix);
	```
=== "Julia"
	```julia
	Op(type::String)
	Op(type::String, site::Int64)
	Op(type::String, sites::Vector{Int64})
	Op(type::String, site::Int64, matrix::Matrix{Float64})
	Op(type::String, sites::Vector{Int64}, matrix::Matrix{Float64})
	Op(type::String, site::Int64, matrix::Matrix{ComplexF64})
	Op(type::String, sites::Vector{Int64}, matrix::Matrix{ComplexF64})
	```

| Parameter | Description                                                   |          |
|:----------|:--------------------------------------------------------------|----------|
| type      | a string which denotes what kind of operator is represented   |          |
| sites     | defines on which site(s) of the lattice the operator acts on. | optional |
| matrix    | defines a matrix which may be needed to describe an operator. | optional |

!!! warning "1-indexing in Julia / 0-indexing in C++"

	To enumerate the sites of an Op, we start counting at 1 in Julia and 0 in C++.

--- 

## Methods

#### isreal

Returns whether an [Op](op.md) is a real operator.

=== "C++"
	```c++
    bool isreal(Op const &op);
	```

=== "Julia"
	```julia
    isreal(op::Op)::Bool
	```
---

#### isapprox

Returns whether two Ops are approximately equal.
	
=== "C++"	
	```c++
	bool isapprox(Op const &op1, OpSum const &op2, double rtol = 1e-12,
	              double atol = 1e-12);
	```

=== "Julia"
	```julia
	isapprox(op1::Op, op2::Op, rtol::Float64=1e-12, atol::Float64=1e-12)::Bool
	```
---

#### to_string (operator<<)

Converts the Op to a readable string representation.
	
=== "C++"	
	```c++
	std::string to_string(Op const &op);
	std::ostream &operator<<(std::ostream &out, Op const &op);
	```

=== "Julia"
	```julia
	to_string(op::Op)::String
	```
---
## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:Op"
	```

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:Op"
	```
