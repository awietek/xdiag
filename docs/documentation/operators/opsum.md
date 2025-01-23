---
title: OpSum
---

Object representing a generic many-body operator by a sum of operators of the form 

$$ \mathcal{O} = \sum_i c_i \mathcal{O}_i. $$

Hence, an OpSum is consists of a sum of pairs given by

1. A coupling constant $c_i$ which is defined by a [Coupling](coupling.md) object. The coupling can either be a string name or a real/complex number.

2. An operator $\mathcal{O}_i$ defined by an [Op](op.md) object.


Generically, an OpSum can thus have coupling constants defined by either strings or numerical real/complex numbers. We call an OpSum **plain** if its couplings are only numerical numbers, and not strings. String couplings can be defined by using the access `operator[]`. If all string coupling constants are defined, the OpSum can be converted to a plain OpSum using the `plain` method shown below.

Thus, OpSums can be defined independently of the numerical values of the coupling constants, e.g. in an input file. Upon execution of the code, these constants can then be set. Most operations in XDiag require the OpSum to be convertible to a plain OpSum.

OpSums can be added and subtracted, as well as multiplied with and divided by a [Scalar](scalar.md). Hence, OpSums carry the mathematical structure of a vector space.

**Sources** [opsum.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/operators/opsum.hpp), [opsum.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/operators/opsum.cpp)

## Constructors

The following constructors create an OpSum with a single pair of coupling and operator. Additional terms can be added using the `+` and `+=` operators explained further below. If no [Coupling](coupling.md) is given, a numerical coefficient of `1.0` is assumed.


=== "C++"	
	```c++
	OpSum(Op const &op);
	OpSum(Coupling const &cpl, Op const &op);
	```

=== "Julia"
	```julia
	OpSum(op::Op)
	OpSum(cpl::Coupling, op::Op)
	```
	
| Parameter | Description                                                                 | Default |
|:----------|:----------------------------------------------------------------------------|---------|
| cpl       | A [Coupling](coupling.md) which is either a string or a real/complex number | 1.0     |
| op        | An [Op](coupling.md) which describes the type of operator                   |         |

Alternatively, an OpSum can also be constructed via the `* operator`, for example:

=== "C++"
	```c++
	auto ops = OpSum();
	for (int i = 0; i<N; ++i) {
      ops += J * Op("SzSz", {i, (i + 1) % N});
    }
	```
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

Creates an OpSum with a single pair of [Coupling](coupling.md) and [Op](op.md) objects.

=== "C++"
	```c++
	OpSum operator*(std::string cpl, Op const &op);
	OpSum operator*(double cpl, Op const &op);
	OpSum operator*(complex cpl, Op const &op);
	OpSum operator*(Scalar cpl, Op const &op);
	OpSum operator*(Coupling const& cpl, Op const &op);
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
	Base.:+(ops::OpSum, ops2::OpSum)
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

#### operator* , operator/ (Scalar muliplication/division)

Multiplies an OpSum $\mathcal{A} = \sum_i a_i \mathcal{A}_i$ with a scalar $b$ to form

$$\mathcal{B} = b \sum_i a_i \mathcal{A}_i$$

=== "C++"
	```c++
	OpSum &operator*=(Scalar const &cpl);
	OpSum &operator/=(Scalar const &cpl);
	OpSum operator*(Scalar const &cpl, OpSum const &op);
	OpSum operator*(OpSum const &op, Scalar const &cpl);
	OpSum operator/(OpSum const &op, Scalar const &cpl);
	```
---

#### operator[]

Sets a coupling constant defined as a string to a numerical value.

=== "C++"
	```c++
	Scalar &operator[](std::string name);
	Scalar const &operator[](std::string name) const;
	```
---

#### constants

Returns a vector of strings with the coupling constants defined, i.e. the strings that define some of the [Coupling](coupling.md) objects.

=== "C++"
	```c++
	std::vector<std::string> constants(OpSum const &ops) const;
	```
---

## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:OpSum"
	```
