---
title: OpSum
---

Object representing a generic many-body operator by a sum of operators of the form 
$$ \mathcal{O} = \sum_i c_i \mathcal{O}_i. $$

**Sources:** [opsum.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/operators/opsum.hpp) · [opsum.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/operators/opsum.cpp)

An OpSum is consists of a sum of pairs given by

1. A coupling constant $c_i$ which is given by a either a string name or a real/complex number.

2. An operator $\mathcal{O}_i$ defined by an [Op](op.md) object.

Generically, an OpSum can thus have coupling constants defined by either strings or numerical real/complex numbers. We call an OpSum **plain** if its couplings are only numerical numbers, and not strings. String couplings can be defined by using the access `operator[]`. If all string coupling constants are defined, the OpSum can be converted to a plain OpSum using the `plain` method shown below.

Thus, OpSums can be defined independently of the numerical values of the coupling constants, e.g. in an input file. Upon execution of the code, these constants can then be set. Most operations in XDiag require the OpSum to be convertible to a plain OpSum.

OpSums can be added and subtracted, as well as multiplied with and divided by a scalar value, i.e. a real or complex number. In addition, two OpSums can be **multiplied** with one another, forming the (generally non-commutative) product of the two operators. Together, these operations turn the set of OpSums into an **algebra**. Complemented by the Hermitian conjugation [hc](hc.md), which is an involution satisfying $(\mathcal{A}\mathcal{B})^\dagger = \mathcal{B}^\dagger\mathcal{A}^\dagger$, the OpSums form an **involutive algebra** ($*$-algebra).

## Constructors

The following constructors create an OpSum with a single pair of coupling and operator. Additional terms can be added using the `+` and `+=` operators explained further below. If no coupling is given, a numerical coefficient of `1.0` is assumed.

=== "Julia"
	```julia
	OpSum(op::Op)
	OpSum(coupling::Float64, op::Op)
	OpSum(coupling::ComplexF64, op::Op)
	OpSum(coupling::String, op::Op)
	```
=== "C++"	
	```c++
	OpSum(Op const &op);
	OpSum(double coupling, Op const &op);
	OpSum(complex coupling, Op const &op);
	OpSum(std::string coupling, Op const &op);
	```
	
| Parameter | Description                                                  | Default |
|:----------|:-------------------------------------------------------------|---------|
| coupling  | A coupling which is either a string or a real/complex number | 1.0     |
| op        | An [Op](op.md) which describes the type of operator          |         |

Alternatively, an OpSum can also be constructed via the `* operator`, for example:

=== "Julia"
	```julia
	ops = OpSum();
	for i in 1:n
		ops += "J" * Op("SzSz", [i, mod1(i+1, N)]);
	```
=== "C++"
	```c++
	auto ops = OpSum();
	for (int i = 0; i<N; ++i) {
		ops += "J" * Op("SzSz", {i, (i + 1) % N});
    }
	```
	
## Complex couplings

XDiag allows all couplings to be complex. A complex coupling always acts as a plain prefactor $c$ multiplying the operator. Under Hermitian conjugation with [hc](hc.md), the coupling is complex-conjugated to $c^*$ (while the operator type is mapped to its Hermitian-conjugate partner, e.g. $\mathrm{hc}(S^+) = S^-$),

$$ \mathrm{hc}(c\,\mathcal{O}) = c^*\, \mathrm{hc}(\mathcal{O}). $$

Consequently, a complex coupling generally turns a Hermitian operator into a non-Hermitian one. This holds for **all** operator types, including hopping (`Hop`, `Hopup`, `Hopdn`) and exchange (`Exchange`) terms: a complex coupling on such a term is *not* invariant under [hc](hc.md). To build a Hermitian operator from a complex coupling, add the Hermitian conjugate explicitly, e.g. `ops + hc(ops)`.

!!! warning "Changed in version 0.5"

	Prior to version 0.5, the hopping (`Hop`, `Hopup`, `Hopdn`) and exchange (`Exchange`) operators were defined to be Hermitian even in the presence of a complex coupling, i.e. the coupling entered as $t\,c^\dagger_i c_j + t^* c^\dagger_j c_i$. This convention was changed, because with it the [OpSum](opsum.md) objects would not form an [algebra](#operator-algebra-product). Since then, a complex coupling is a plain prefactor as described above. To obtain the antisymmetric (non-Hermitian) hopping and exchange terms, the dedicated operator types [`HopAsym`](operator_types.md) and [`ExchangeAsym`](operator_types.md) (as well as the spin-resolved `HopupAsym` and `HopdnAsym`) were introduced.


## Methods

#### plain

Converts an OpSum with possible string couplings to an OpSum with purely numerical real/complex couplings.

=== "Julia"
	```julia
	plain(ops::OpSum)::OpSum
	```
=== "C++"
	```c++
	OpSum plain(OpSum const &ops) const;
	```

#### operator* (Creation)

Creates an OpSum with a single pair of coupling constant and an [Op](op.md) object.

=== "Julia"
	```julia
	Base.:*(coupling::Int64, op::Op)::OpSum
	Base.:*(coupling::Float64, op::Op)::OpSum
	Base.:*(coupling::ComplexF64, op::Op)::OpSum
	Base.:*(coupling::String, op::Op)::OpSum
	```
=== "C++"
	```c++
	OpSum operator*(int64_t coupling, Op const &op);
	OpSum operator*(double coupling, Op const &op);
	OpSum operator*(complex coupling, Op const &op);
	OpSum operator*(std::string coupling, Op const &op);
	```

#### operator+ / operator +=

Adds two OpSum objects $\mathcal{A} = \sum_i a_i \mathcal{A}_i$ and $\mathcal{B} = \sum_i b_i \mathcal{B}_i$ to for the sum of the two operators,
	$$ \mathcal{A} + \mathcal{B} = \sum_i a_i \mathcal{A}_i + \sum_i b_i \mathcal{B}_i$$

=== "Julia"
	```julia
	Base.:+(ops1::OpSum, ops2::OpSum)::OpSum
	```
=== "C++"
	```c++
	OpSum &operator+=(OpSum const &ops);
	OpSum operator+(OpSum const &ops) const;
	```
	
#### operator- / operator -=

Subtracts to OpSum objects.

=== "Julia"
	```julia
	Base.:-(ops::OpSum, ops2::OpSum)::OpSum
	```
=== "C++"
	```c++
	OpSum &operator-=(OpSum const &ops);
	OpSum operator-(OpSum const &ops) const;
	```

#### operator* , operator/ (scalar muliplication/division)

Multiplies an OpSum $\mathcal{A} = \sum_i a_i \mathcal{A}_i$ with a scalar $b$ to form

$$\mathcal{B} = b \sum_i a_i \mathcal{A}_i$$

=== "Julia"
	```julia
	Base.:*(coupling::Int64, ops::OpSum)::OpSum
	Base.:*(coupling::Float64, ops::OpSum)::OpSum
	Base.:*(coupling::ComplexF64, ops::OpSum)::OpSum
	Base.:*(ops::OpSum, coupling::Int64)::OpSum
	Base.:*(ops::OpSum, coupling::Float64)::OpSum
	Base.:*(ops::OpSum, coupling::ComplexF64)::OpSum
	Base.:/(ops::OpSum, coupling::Int6464)::OpSum
	Base.:/(ops::OpSum, coupling::Float64)::OpSum
	Base.:/(ops::OpSum, coupling::ComplexF64)::OpSum
	```
=== "C++"
	```c++
	OpSum &operator*=(int64_t scalar);
	OpSum &operator*=(double scalar);
    OpSum &operator*=(complex scalar);
	OpSum &operator/=(int64_t scalar);
	OpSum &operator/=(double scalar);
	OpSum &operator/=(complex scalar);
  
	OpSum operator*(int64_t scalar, OpSum const &op);
	OpSum operator*(double scalar, OpSum const &op);
	OpSum operator*(complex scalar, OpSum const &op);
	OpSum operator*(OpSum const &op, double scalar);
	OpSum operator*(OpSum const &op, int64_t scalar);
	OpSum operator*(OpSum const &op, complex scalar);
	OpSum operator/(OpSum const &op, int64_t scalar);
	OpSum operator/(OpSum const &op, double scalar);
	OpSum operator/(OpSum const &op, complex scalar);
	```

#### operator* (algebra product)

Multiplies two OpSum objects $\mathcal{A} = \sum_i a_i \mathcal{A}_i$ and $\mathcal{B} = \sum_j b_j \mathcal{B}_j$ to form their operator product, distributing over the sums,

$$ \mathcal{A}\,\mathcal{B} = \sum_{i,j} (a_i b_j)\, \mathcal{A}_i \mathcal{B}_j. $$

The product is generally **non-commutative**, reflecting that quantum mechanical operators need not commute. It can be used to build composite operators, e.g. correlation functions or powers of a Hamiltonian, from elementary building blocks.

=== "Julia"
	```julia
	Base.:*(ops1::OpSum, ops2::OpSum)::OpSum
	```
=== "C++"
	```c++
	OpSum &operator*=(OpSum const &ops);
	OpSum operator*(OpSum const &ops) const;
	```

#### operator[]

Sets a coupling constant defined as a string to a numerical value.

=== "Julia"
	```julia
	Base.setindex!(ops::OpSum, cpl::Int64, name::String)
	Base.setindex!(ops::OpSum, cpl::Float64, name::String)
	Base.setindex!(ops::OpSum, cpl::ComplexF64, name::String)
	```
=== "C++"
	```c++
	Scalar &operator[](std::string name);
	```

#### isreal

Returns whether an [OpSum](opsum.md) is a real operator.

=== "Julia"
	```julia
    isreal(ops::OpSum)::Bool
	```
=== "C++"
	```c++
    bool isreal(OpSum const &ops);
	```

#### to_string (operator<<)

Converts the OpSum to a readable string representation.

=== "Julia"
	```julia
	to_string(ops::OpSum)::String
	```
=== "C++"	
	```c++
	std::string to_string(OpSum const &ops);
	std::ostream &operator<<(std::ostream &out, OpSum const &ops);
	```


	
## Usage Example

=== "Julia"
	```julia
	--8<-- "examples/usage_examples/main.cpp:opsum"
	```
=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:OpSum"
	```
