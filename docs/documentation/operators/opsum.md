---
title: OpSum
---

A sum of operators of the form $$ \mathcal{O} = \sum_i c_i \mathcal{O}_i. $$
The coupling constants $c_i$ are given by a [Coupling](coupling.md) object, which can either be
a string or a real/complex number. The operators $\mathcal{O}_i$ are defined by an [Op](coupling.md) object.

**Source** [opsum.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/operators/opsum.hpp)

## Constructors

=== "C++"	
	```c++
	OpSum() = default;
	OpSum(Op const &op);
	OpSum(Coupling const &cpl, Op const &op);
	```

| Parameter | Description                                                                 | Default |
|:----------|:----------------------------------------------------------------------------|---------|
| cpl       | A [Coupling](coupling.md) which is either a string or a real/complex number | 1.0     |
| Op        | An [Op](coupling.md) which describes the type of operator                   |         |

---

## Methods

#### operator+

Adds two OpSum objects $\mathcal{A} = \sum_i a_i \mathcal{A}_i$ and $\mathcal{B} = \sum_i b_i \mathcal{B}_i$ to for the sum of the two operators,
	$$ \mathcal{A} + \mathcal{B} = \sum_i a_i \mathcal{A}_i + \sum_i b_i \mathcal{B}_i$$

=== "C++"
	```c++
	void operator+=(OpSum const &ops);
	OpSum operator+(OpSum const &ops) const;
	```
---

#### operator*

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

## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:opsum"
	```
