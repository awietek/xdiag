---
title: OpSum
---

Object representing a generic many-body operator by a sum of operators of the form 

$$ \mathcal{O} = \sum_i c_i \mathcal{O}_i. $$

**Sources**<br>
[opsum.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/operators/opsum.hpp)<br>
[opsum.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/operators/opsum.cpp)<br>
[opsum.jl](https://github.com/awietek/XDiag.jl/blob/main/src/operators/opsum.jl)

--- 

An OpSum is consists of a sum of pairs given by

1. A coupling constant $c_i$ which is given by a either a string name or a real/complex number.

2. An operator $\mathcal{O}_i$ defined by an [Op](op.md) object.

Generically, an OpSum can thus have coupling constants defined by either strings or numerical real/complex numbers. We call an OpSum **plain** if its couplings are only numerical numbers, and not strings. String couplings can be defined by using the access `operator[]`. If all string coupling constants are defined, the OpSum can be converted to a plain OpSum using the `plain` method shown below.

Thus, OpSums can be defined independently of the numerical values of the coupling constants, e.g. in an input file. Upon execution of the code, these constants can then be set. Most operations in XDiag require the OpSum to be convertible to a plain OpSum.

OpSums can be added and subtracted, as well as multiplied with and divided by a scalar value, i.e. a real or complex number. Hence, OpSums carry the mathematical structure of a vector space.

---

## Constructors

The following constructors create an OpSum with a single pair of coupling and operator. Additional terms can be added using the `+` and `+=` operators explained further below. If no coupling is given, a numerical coefficient of `1.0` is assumed.


=== "C++"	
	```c++
	OpSum(Op const &op);
	OpSum(double coupling, Op const &op);
	OpSum(complex coupling, Op const &op);
	OpSum(std::string coupling, Op const &op);
	```

=== "Julia"
	```julia
	OpSum(op::Op)
	OpSum(coupling::Float64, op::Op)
	OpSum(coupling::ComplexF64, op::Op)
	OpSum(coupling::String, op::Op)
	```
	
| Parameter | Description                                                  | Default |
|:----------|:-------------------------------------------------------------|---------|
| coupling  | A coupling which is either a string or a real/complex number | 1.0     |
| op        | An [Op](op.md) which describes the type of operator          |         |

Alternatively, an OpSum can also be constructed via the `* operator`, for example:

=== "C++"
	```c++
	auto ops = OpSum();
	for (int i = 0; i<N; ++i) {
		ops += "J" * Op("SzSz", {i, (i + 1) % N});
    }
	```
=== "Julia"
	```julia
	ops = OpSum();
	for i in 1:n
		ops += "J" * Op("SzSz", [i, mod1(i+1, N)]);
	```
	
---

## Complex couplings

XDiag allows all couplings to be complex. Depending on the [operator type](operator_types.md) 
a complex coupling can have two meanings:

1. A complex prefactor $c$ which upon hermitian conjugation with [hc](hc.md) gets 
   conjugated to $c^\star$. This is the case for the following interaction types:   
   `HubbardU`, `Cdagup`, `Cdagdn`, `Cup`, `Cdn`, `Nup`, `Ndn`, `Ntot`, `NtotNtot`, 
   `SdotS`, `SzSz`, `Sz`, `S+`, `S-`, `ScalarChirality`, `tJSzSz`, `tJSdotS`, 
   `Matrix`  
   Thus, a complex coupling can turn a Hermitian operator to a non-Hermitian operator.

2. The coupling is part of the definition of the operator. For, example a hopping 
   operator of the form   
   $$ ( t c^\dagger_{i\sigma}c_{j\sigma} + \textrm{h.c.})  = ( t c^\dagger_{i\sigma}c_{j\sigma} + t^\star c^\dagger_{j\sigma}c_{i\sigma}) $$  
   A complex coupling $t$ gives the hopping a phase, but the overall operator remains
   Hermitian and, thus, invariant under [hc](hc.md). This holds for the types 
   `Hop`, `Hopup`, `Hopdn`, `Exchange`. In the latter case, complex spin exchange 
   `Exchange` is defined as,
   $$ \frac{1}{2}( J S^+_i S^-_j + J^\star S^-_iS^+_j)$$
	
---

## Methods

#### plain

Converts an OpSum with possible string couplings to an OpSum with purely numerical real/complex couplings.

=== "C++"
	```c++
	OpSum plain(OpSum const &ops) const;
	```
	
=== "Julia"
	```julia
	plain(ops::OpSum)
	```
---

#### operator* (Creation)

Creates an OpSum with a single pair of coupling constant and an [Op](op.md) object.

=== "C++"
	```c++
	OpSum operator*(double coupling, Op const &op);
	OpSum operator*(complex coupling, Op const &op);
	OpSum operator*(std::string coupling, Op const &op);
	```
=== "Julia"
	```julia
	Base.:*(coupling::Float64, op::Op)
	Base.:*(coupling::ComplexF64, op::Op)
	Base.:*(coupling::String, op::Op)
	```
---

#### operator+ / operator +=

Adds two OpSum objects $\mathcal{A} = \sum_i a_i \mathcal{A}_i$ and $\mathcal{B} = \sum_i b_i \mathcal{B}_i$ to for the sum of the two operators,
	$$ \mathcal{A} + \mathcal{B} = \sum_i a_i \mathcal{A}_i + \sum_i b_i \mathcal{B}_i$$

=== "C++"
	```c++
	OpSum &operator+=(OpSum const &ops);
	OpSum operator+(OpSum const &ops) const;
	```
	
=== "Julia"
	```julia
	Base.:+(ops1::OpSum, ops2::OpSum)
	```

---
#### operator- / operator -=

Subtracts to OpSum objects.

=== "C++"
	```c++
	OpSum &operator-=(OpSum const &ops);
	OpSum operator-(OpSum const &ops) const;
	```
	
=== "Julia"
	```julia
	Base.:-(ops::OpSum, ops2::OpSum)
	```
---

#### operator* , operator/ (scalar muliplication/division)

Multiplies an OpSum $\mathcal{A} = \sum_i a_i \mathcal{A}_i$ with a scalar $b$ to form

$$\mathcal{B} = b \sum_i a_i \mathcal{A}_i$$

=== "C++"
	```c++
	OpSum &operator*=(double scalar);
    OpSum &operator*=(complex scalar);
	OpSum &operator/=(double scalar);
	OpSum &operator/=(complex scalar);
  
	OpSum operator*(double scalar, OpSum const &op);
	OpSum operator*(complex scalar, OpSum const &op);
	OpSum operator*(OpSum const &op, double scalar);
	OpSum operator*(OpSum const &op, complex scalar);
	OpSum operator/(OpSum const &op, double scalar);
	OpSum operator/(OpSum const &op, complex scalar);
	```

=== "Julia"
	```julia
	Base.:*(coupling::Float64, ops::OpSum)
	Base.:*(coupling::ComplexF64, ops::OpSum)
	Base.:*(ops::OpSum, coupling::Float64)
	Base.:*(ops::OpSum, coupling::ComplexF64)
	Base.:/(ops::OpSum, coupling::Float64)
	Base.:/(ops::OpSum, coupling::ComplexF64)
	```

---

#### operator[]

Sets a coupling constant defined as a string to a numerical value.

=== "C++"
	```c++
	Scalar &operator[](std::string name);
	```
=== "Julia"
	```julia
	Base.setindex!(ops::OpSum, cpl::Float64, name::String)
	Base.setindex!(ops::OpSum, cpl::ComplexF64, name::String)
	```	
---

#### constants

Returns a vector of strings with the coupling constants defined, i.e. the strings that define some of the coupling constants.

=== "C++"
	```c++
	std::vector<std::string> constants(OpSum const &ops) const;
	```
=== "Julia"
	```c++
	constants(ops::OpSum)
	```
	
---

#### isreal

Returns whether an [OpSum](opsum.md) is a real operator.

=== "C++"
	```c++
    bool isreal(OpSum const &ops);
	```

=== "Julia"
	```julia
    isreal(ops::OpSum)
	```
---

#### isapprox

Returns whether two OpSums are approximately equal.
	
=== "C++"	
	```c++
	bool isapprox(OpSum const &ops1, OpSum const &ops2, double rtol = 1e-12,
	              double atol = 1e-12);
	```

=== "Julia"
	```julia
	isapprox(ops1::OpSum, ops2::OpSum, rtol::Float64=1e-12, atol::Float64=1e-12)
	```
---

#### to_string (operator<<)

Converts the OpSum to a readable string representation.
	
=== "C++"	
	```c++
	std::string to_string(OpSum const &ops);
	std::ostream &operator<<(std::ostream &out, OpSum const &ops);
	```

=== "Julia"
	```julia
	to_string(ops::OpSum)
	```
---
	
## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:OpSum"
	```
