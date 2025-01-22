---
title: tJ
---

Representation of a block in a  $t-J$ type Hilbert space. 

**Source** [tj.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/blocks/tj.hpp)

## Constructors

=== "Julia"
	```julia
	tJ(nsites::Integer, nup::Integer, ndn::Integer)
	tJ(nsites::Integer, nup::Integer, ndn::Integer, 
	   group::PermutationGroup, irrep::Representation)
	```

=== "C++"	
	```c++
    tJ(int64_t nsites, int64_t nup, int64_t ndn);
    tJ(int64_t nsites, int64_t nup, int64_t ndn, 
	   PermutationGroup group, Representation irrep);
	```


| Name    | Description                                                                                |   |
|:--------|:-------------------------------------------------------------------------------------------|---|
| nsites | number of sites (integer)                                                                  |   |
| nup    | number of "up" electrons (integer)                                                         |   |
| ndn    | number of "dn" electrons (integer)                                                         |   |
| group   | [PermutationGroup](../symmetries/permutation_group.md) defining the permutation symmetries |   |
| irrep   | Irreducible [Representation](../symmetries/representation.md)  of the symmetry group       |   |

## Iteration

An tJ block can be iterated over, where at each iteration a [ProductState](../states/product_state.md) representing the corresponding basis state is returned.

=== "Julia"
	```julia
	block = tJ(4, 2, 1)
	for pstate in block
		@show pstate, index(block, pstate) 
	end
	```

=== "C++"	
	```c++
    auto block = tJ(4, 2, 1);
	for (auto pstate : block) {
		Log("{} {}", to_string(pstate), block.index(pstate));
	}
	```

## Methods

!!! method "index"

	Returns the index of a given [ProductState](../states/product_state.md) in the basis of the tJ block.

	=== "Julia"
		```julia
		index(block::tJ, pstate::ProductState)
		```

	=== "C++"	
		```c++
		int64_t index(ProductState const &pstate) const;
		```
	!!! warning "1-indexing"
		In the C++ version, the index count starts from "0" whereas in Julia the index count starts from "1".



!!! method "nsites"

	Returns the number of sites of the block.

	=== "Julia"
		```julia
		nsites(block::tJ)
		```

	=== "C++"	
		```c++
		int64_t nsites() const;
		```

!!! method "nup"

	Returns the number of "up" electrons.

	=== "Julia"
		```julia
		nup(block::tJ)
		```

	=== "C++"	
		```c++
		int64_t nup() const;
		```


!!! method "ndn"

	Returns the number of "down" electrons.

	=== "Julia"
		```julia
		ndn(block::tJ)
		```

	=== "C++"	
		```c++
		int64_t ndn() const;
		```

!!! method "permutation_group"

	Returns the [PermutationGroup]("../symmetries/permutation_group.md") of the block, if defined.

	=== "Julia"
		```julia
		permutation_group(block::tJ)
		```

	=== "C++"	
		```c++
	    PermutationGroup permutation_group() const;
		```


!!! method "irrep"

	Returns the [Representation]("../symmetries/representation.md") of the block, if defined.

	=== "Julia"
		```julia
	    irrep(block::tJ)
		```

	=== "C++"	
		```c++
	    Representation irrep() const;
		```


!!! method "size"
	Returns the size of the block, i.e. its dimension.

	=== "Julia"
		```julia
		size(block::tJ)
		```

	=== "C++"	
		```c++
		int64_t size() const;
		```

!!! method "dim"
	Returns the dimension of the block, same as "size" for non-distributed blocks.

	=== "Julia"
		```julia
		dim(block::tJ)
		```

	=== "C++"	
		```c++
		int64_tdim() const;
		```
		
!!! method "isreal"
	Returns whether the block can be used with real arithmetic. 
	Complex arithmetic is needed when a
	[Representation](../symmetries/representation.md) is genuinely complex.

	=== "Julia"
		```julia
	    isreal(block::tJ; precision::Real=1e-12)
		```

	=== "C++"	
		```c++
		int64_t isreal(double precision = 1e-12) const;
		```

## Usage Example

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:tJ"
	```

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:tJ"
	```

