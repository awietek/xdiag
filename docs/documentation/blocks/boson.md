---
title: Boson
---

A block in a Hilbert space of bosons with a fixed local dimension $d$, i.e. each site can be occupied by $0, 1, \ldots, d-1$ particles.

Since a $d$-level local degree of freedom equally describes a spin $S = (d-1)/2$, the same block also represents general spin-$S$ systems. For this reason `Boson` is also available under the alias `Spin`, and both the bosonic ladder operators and the spin-$S$ operators are defined on it (see [Operators](#operators) below).

**Sources:** [boson.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/blocks/boson.hpp) · [boson.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/blocks/boson.cpp)

## Constructors

=== "C++"
	```c++
	Boson(int64_t nsites, int64_t d);
	Boson(int64_t nsites, int64_t d, int64_t number);
	Boson(int64_t nsites, int64_t d, Representation const &irrep);
	Boson(int64_t nsites, int64_t d, int64_t number, Representation const &irrep);
	```
=== "Julia"
	```julia
	Boson(nsites::Int64, d::Int64)
	Boson(nsites::Int64, d::Int64, number::Int64)
	Boson(nsites::Int64, d::Int64, irrep::Representation)
	Boson(nsites::Int64, d::Int64, number::Int64, irrep::Representation)
	```

| Name    | Description                                                                          | Default |
|:--------|:-------------------------------------------------------------------------------------|---------|
| nsites  | number of sites (integer)                                                            |         |
| d       | local dimension, i.e. number of local states $0,\ldots,d-1$ ($2 \le d \le 256$)      |         |
| number  | total number of bosons (integer)                                                     |         |
| irrep   | Irreducible [Representation](../symmetries/representation.md)  of the symmetry group |         |

If the `number` argument is omitted, the block contains all boson numbers. The spin-$S$ degree of freedom corresponds to a local dimension $d = 2S + 1$; the `Spin` alias can be used interchangeably with `Boson`.

## Local configurations

Each site of a Boson block carries a local dimension $d$. In a [ProductState](../states/product_state.md), the local configuration of a site is simply the **occupation number**, i.e. an integer in $0, 1, \ldots, d-1$. In the spin-$S$ interpretation, an occupation $n$ corresponds to the magnetic quantum number $m = n - S$ with $S = (d-1)/2$.

## Operators

Two families of operators are defined on a Boson block. The bosonic ladder operators act on the occupation numbers, with $a_i |n\rangle = \sqrt{n}\,|n-1\rangle$, $a^\dagger_i |n\rangle = \sqrt{n+1}\,|n+1\rangle$ (truncated at $n = d-1$), and $n_i = a^\dagger_i a_i$.

| Type       | Description                        | Formula                                            | No. of sites |
|:-----------|:-----------------------------------|:---------------------------------------------------|:------------:|
| `Adag`     | bosonic creation operator          | $a^\dagger_i$                                      | 1            |
| `A`        | bosonic annihilation operator      | $a_i$                                              | 1            |
| `N`        | number operator                    | $n_i = a^\dagger_i a_i$                            | 1            |
| `Hop`      | hopping term                       | $-(a^\dagger_i a_j + a^\dagger_j a_i)$             | 2            |
| `HopAsym`  | antisymmetric hopping term         | $-(a^\dagger_i a_j - a^\dagger_j a_i)$             | 2            |
| `HubbardU` | on-site interaction                | $\frac{1}{2}\sum_i n_i (n_i - 1)$                  | 0            |
| `TotalN`   | total number of bosons             | $\sum_i n_i$                                        | 0            |

The spin-$S$ operators (with $S = (d-1)/2$) act on the same local states. Here $\mathbf{S}_i = (S^x_i, S^y_i, S^z_i)$ and $S^\pm_i = S^x_i \pm i S^y_i$.

| Type              | Description                        | Formula                                                                       | No. of sites |
|:------------------|:-----------------------------------|:------------------------------------------------------------------------------|:------------:|
| `Sz`              | local magnetic moment ($z$)        | $S^z_i$                                                                       | 1            |
| `Sx`              | local magnetic moment ($x$)        | $S^x_i$                                                                       | 1            |
| `Sy`              | local magnetic moment ($y$)        | $S^y_i$                                                                       | 1            |
| `S+`              | spin raising operator              | $S^+_i$                                                                       | 1            |
| `S-`              | spin lowering operator             | $S^-_i$                                                                       | 1            |
| `SdotS`           | Heisenberg interaction             | $\mathbf{S}_i \cdot \mathbf{S}_j$                                             | 2            |
| `SzSz`            | Ising interaction                  | $S^z_i S^z_j$                                                                 | 2            |
| `Exchange`        | spin exchange interaction          | $\frac{1}{2}(S^+_i S^-_j + S^-_i S^+_j)$                                      | 2            |
| `ExchangeAsym`    | antisymmetric exchange             | $\frac{1}{2}(S^+_i S^-_j - S^-_i S^+_j)$                                      | 2            |
| `ScalarChirality` | scalar chirality                   | $\mathbf{S}_i \cdot (\mathbf{S}_j \times \mathbf{S}_k)$                       | 3            |
| `Matrix`          | generic operator via a matrix      | user-defined matrix on the $d^n$-dimensional local space of $n$ sites          | any          |
| `Id`              | identity                           | $\mathbb{1}$                                                                  | 0            |

For a full description of all operator types, see the [operator types](../operators/operator_types.md) page.

## Iteration

A Boson block can be iterated over, where at each iteration a [ProductState](../states/product_state.md) representing the corresponding basis state is returned.

=== "C++"
	```c++
	auto block = Boson(4, 3, 2);   // 4 sites, local dimension 3, 2 bosons
	for (auto pstate : block) {
	  Log("{} {}", to_string(pstate, block), block.index(pstate));
	}
	```
=== "Julia"
	```julia
	block = Boson(4, 3, 2)   # 4 sites, local dimension 3, 2 bosons
	for pstate in block
	    @show to_string(pstate, block), index(block, pstate)
	end
	```

## Methods

#### index
Returns the index of a given [ProductState](../states/product_state.md) in the basis of the Boson block.

=== "C++"
	```c++
	int64_t index(Boson const &block, ProductState const &pstate);
	```
=== "Julia"
	```julia
	index(block::Boson, pstate::ProductState)::Int64
	```

!!! warning "1-indexing"
	In the C++ version, the index count starts from "0" whereas in Julia the index count starts from "1".

#### nsites
Returns the number of sites of the block.

=== "C++"
	```c++
	int64_t nsites(Boson const &block);
	```
=== "Julia"
	```julia
	nsites(block::Boson)::Int64
	```

#### d
Returns the local dimension of the block.

=== "C++"
	```c++
	int64_t d(Boson const &block);
	```
=== "Julia"
	```julia
	d(block::Boson)::Int64
	```

#### size
Returns the size of the block, i.e. its dimension.

=== "C++"
	```c++
	int64_t size(Boson const &block);
	```
=== "Julia"
	```julia
	size(block::Boson)::Int64
	```

#### dim
Returns the dimension of the block, same as "size" for non-distributed blocks.

=== "C++"
	```c++
	int64_t dim(Boson const &block);
	```
=== "Julia"
	```julia
	dim(block::Boson)::Int64
	```

#### isreal
Returns whether the block can be used with real arithmetic. Complex arithmetic is needed when a [Representation](../symmetries/representation.md) is genuinely complex.

=== "C++"
	```c++
	bool isreal(Boson const &block);
	```
=== "Julia"
	```julia
	isreal(block::Boson)::Bool
	```

## Usage Example

A Bose-Hubbard chain with nearest-neighbor hopping and on-site interaction.

=== "Julia"
	```julia
	--8<-- "examples/usage_examples/main.jl:Boson"
	```
=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:Boson"
	```
