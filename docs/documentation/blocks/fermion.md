---
title: Fermion
---

A block in a Hilbert space of spinless fermions.

**Sources:** [fermion.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/blocks/fermion.hpp) · [fermion.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/blocks/fermion.cpp)

## Constructors

=== "Julia"
	```julia
	Fermion(nsites::Int64)
	Fermion(nsites::Int64, number::Int64)
	Fermion(nsites::Int64, irrep::Representation)
	Fermion(nsites::Int64, number::Int64, irrep::Representation)
	```
=== "C++"
	```c++
	Fermion(int64_t nsites);
	Fermion(int64_t nsites, int64_t number);
	Fermion(int64_t nsites, Representation const &irrep);
	Fermion(int64_t nsites, int64_t number, Representation const &irrep);
	```


| Name    | Description                                                                          | Default |
|:--------|:-------------------------------------------------------------------------------------|---------|
| nsites  | number of sites (integer)                                                            |         |
| number  | number of fermions (integer)                                                         |         |
| irrep   | Irreducible [Representation](../symmetries/representation.md)  of the symmetry group |         |

If the `number` argument is omitted, the block contains all fermion numbers from $0$ to `nsites`.

## Local configurations

Each site of a Fermion block carries a local dimension $d=2$. In a [ProductState](../states/product_state.md), the local configuration of every site is given by an integer with the following meaning:

| Integer | Configuration     | Symbol |
|:--------|:------------------|:-------|
| `0`     | empty             | ○      |
| `1`     | occupied fermion  | ●      |

## Operators

The following operator types can be used on a Fermion block. Here $c^\dagger_i$ and $c_i$ denote the fermionic creation and annihilation operators on site $i$, and $n_i = c^\dagger_i c_i$ is the number operator.

| Type      | Description                 | Formula                                | No. of sites |
|:----------|:----------------------------|:---------------------------------------|:------------:|
| `Cdag`    | creation operator           | $c^\dagger_i$                          | 1            |
| `C`       | annihilation operator       | $c_i$                                  | 1            |
| `N`       | number operator             | $n_i = c^\dagger_i c_i$                | 1            |
| `NN`      | density-density interaction | $n_i n_j$                              | 2            |
| `Hop`     | hopping term                | $-(c^\dagger_i c_j + c^\dagger_j c_i)$ | 2            |
| `HopAsym` | antisymmetric hopping term  | $-(c^\dagger_i c_j - c^\dagger_j c_i)$ | 2            |
| `TotalN`  | total number of fermions    | $\sum_i n_i$                           | 0            |
| `Id`      | identity                    | $\mathbb{1}$                           | 0            |

For a full description of all operator types, see the [operator types](../operators/operator_types.md) page. The Jordan-Wigner sign convention used for the fermionic operators is described in the [Hilbert spaces](../../user_guide/03-hilbert-spaces.md#normal-ordering-of-fermionic-blocks) section of the user guide.

## Iteration

A Fermion block can be iterated over, where at each iteration a [ProductState](../states/product_state.md) representing the corresponding basis state is returned.

=== "Julia"
	```julia
	block = Fermion(4, 2)
	for pstate in block
	    @show to_string(pstate, block), index(block, pstate)
	end
	```
=== "C++"
	```c++
	auto block = Fermion(4, 2);
	for (auto pstate : block) {
	  Log("{} {}", to_string(pstate, block), block.index(pstate));
	}
	```

## Methods

#### index
Returns the index of a given [ProductState](../states/product_state.md) in the basis of the Fermion block.

=== "Julia"
	```julia
	index(block::Fermion, pstate::ProductState)::Int64
	```
=== "C++"
	```c++
	int64_t Fermion::index(ProductState const &pstate);
	```

!!! warning "1-indexing"
	In the C++ version, the index count starts from "0" whereas in Julia the index count starts from "1".

#### nsites
Returns the number of sites of the block.

=== "Julia"
	```julia
	nsites(block::Fermion)::Int64
	```
=== "C++"
	```c++
	int64_t nsites(Fermion const &block);
	```

#### size
Returns the size of the block, i.e. its dimension.

=== "Julia"
	```julia
	size(block::Fermion)::Int64
	```
=== "C++"
	```c++
	int64_t size(Fermion const &block);
	```

#### dim
Returns the dimension of the block, same as "size" for non-distributed blocks.

=== "Julia"
	```julia
	dim(block::Fermion)::Int64
	```
=== "C++"
	```c++
	int64_t dim(Fermion const &block);
	```

#### isreal
Returns whether the block can be used with real arithmetic. Complex arithmetic is needed when a [Representation](../symmetries/representation.md) is genuinely complex.

=== "Julia"
	```julia
	isreal(block::Fermion)::Bool
	```
=== "C++"
	```c++
	bool isreal(Fermion const &block);
	```

## Usage Example

A spinless fermion chain with nearest-neighbor hopping and density-density repulsion.

=== "Julia"
	```julia
	--8<-- "examples/usage_examples/main.jl:Fermion"
	```
=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:Fermion"
	```

