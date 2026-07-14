---
title: States
---

# States

Quantum states $| \psi \rangle$ are represented in XDiag by using a [State](../documentation/states/state.md) object.

## Creating states

The most direct way to create a state is to hand a block together with a vector holding the coefficients on the computational basis of that block. The length of the vector must equal the dimension of the block. The coefficients can be real or complex, given as an `arma::vec` or `arma::cx_vec` in C++, respectively as a `Vector{Float64}` or `Vector{ComplexF64}` in Julia. XDiag automatically creates a real or complex state accordingly.

=== "Julia"
	```julia
	--8<-- "examples/user_guide/main.jl:usage_guide_stat1"
	```
=== "C++"
	```c++
	--8<-- "examples/user_guide/main.cpp:usage_guide_stat1"
	```

A state with zero coefficients is created either implicitly by calling the constructor of `State` with a given block, or explicitly by calling the [zero_state](../documentation/states/create_state.md/#zero_state) function. In both cases, the parameter `real` is optional, can be omitted, and defaults to `true`; it controls whether the state holds real (double precision) or complex (double precision) coefficients.

=== "Julia"
	```julia
	--8<-- "examples/user_guide/main.jl:usage_guide_stat2"
	```
=== "C++"
	```c++
	--8<-- "examples/user_guide/main.cpp:usage_guide_stat2"
	```

We can also create product states. A product state is specified by its per-site local quantum numbers, given as **integers**, following the same convention as the [ProductState](../documentation/states/product_state.md) configurations discussed in the [Hilbert spaces](03-hilbert-spaces.md) section. For a [Spinhalf](../documentation/blocks/spinhalf.md) block, for example, `0` denotes a $\downarrow$-spin and `1` an $\uparrow$-spin. The meaning of the integers for every block type is documented on the respective [block pages](../documentation/blocks/spinhalf.md).

=== "Julia"
	```julia
	--8<-- "examples/user_guide/main.jl:usage_guide_stat3"
	```
=== "C++"
	```c++
	--8<-- "examples/user_guide/main.cpp:usage_guide_stat3"
	```

Finally, a random state with normal $\mathcal{N}(0, 1)$ distributed coefficients is created using the [random_state](../documentation/states/create_state.md/#random_state) function.

=== "Julia"
	```julia
	--8<-- "examples/user_guide/main.jl:usage_guide_stat_rand"
	```
=== "C++"
	```c++
	--8<-- "examples/user_guide/main.cpp:usage_guide_stat_rand"
	```

The random number generation is **deterministic** and controlled by the integer `seed` argument: calling [random_state](../documentation/states/create_state.md/#random_state) with the same `seed` (and the same block) always produces exactly the same state. This is essential for reproducibility. To obtain a different random state, simply choose a different `seed`. The `seed` argument is optional and defaults to a fixed value.

## Norms and overlaps

The $2$-norm $\parallel  | \psi \rangle\parallel_2$ and dot product $\langle \psi_1 | \psi_2 \rangle$ of states can easily be computed using the [norm](../documentation/states/algebra.md/#norm) and [dot/dotC](../documentation/states/algebra.md/#dot) functions.

=== "Julia"
	```julia
	--8<-- "examples/user_guide/main.jl:usage_guide_stat4"
	```
=== "C++"
	```c++
	--8<-- "examples/user_guide/main.cpp:usage_guide_stat4"
	```

The function [dotC](../documentation/states/algebra.md/#dot) is only available in C++, and returns a complex (double precision) number whenever one of the involved states is complex. This is necessary, as the return type of a function must be known at compile time in C++, whereas Julia permits dynamic typing.

## Accessing coefficients

The coefficients of a given state can be retrieved using the [vector/vectorC](../documentation/states/state.md/#vectorvectorc) functions.

=== "Julia"
	```julia
	--8<-- "examples/user_guide/main.jl:usage_guide_stat5"
	```
=== "C++"
	```c++
	--8<-- "examples/user_guide/main.cpp:usage_guide_stat5"
	```

Again, the function [vectorC](../documentation/states/state.md/#vectorvectorc) only exists in the C++ version since the return type needs to be known at compile time. In Julia, the type of the vector is decided at runtime.

## Applying operators

Finally, we can apply an operator [OpSum](../documentation/operators/opsum.md) to a state using the [apply](../documentation/kernels/apply.md) function.

=== "Julia"
	```julia
	--8<-- "examples/user_guide/main.jl:usage_guide_stat6"
	```
=== "C++"
	```c++
	--8<-- "examples/user_guide/main.cpp:usage_guide_stat6"
	```

Importantly, if the block of the [State](../documentation/states/state.md) object has a well-defined quantum number, for example, a conserved particle number, XDiag will automatically detect the quantum number of the resulting state or report an error if the operator does not have a well-defined quantum number. This could be the case, for example, when applying a raising or lowering operator on a particle number conserving state. The [apply](../documentation/kernels/apply.md) function acts on a state without creating a matrix representation of the operator, sometimes referred to as *on-the-fly* matrix application.
