### States

Quantum states $| \psi \rangle$ are represented in XDiag by using a [State](documentation/states/state.md) object. A state with zero coefficients is created either implicitly by calling the constructor of `State` with a given block, or explicitly by calling the [zero_state](documentation/states/create_state.md/#zero_state) function.

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:usage_guide_stat1"
	```
	
=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:usage_guide_stat1"
	```

We hereby create a state with real (double precision) coefficients or complex (double precision) coefficients. The parameter `real` is optional, can be omitted and defaults to `true`. A state with a given vector of coefficients can also be created.

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:usage_guide_stat2"
	```
	
=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:usage_guide_stat2"
	```

Moreover, we can create product states as well as random states (with normal $\mathcal{N}(0, 1)$ distributed coefficients).

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:usage_guide_stat3"
	```
	
=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:usage_guide_stat3"
	```

The $2$-norm $\parallel  | \psi \rangle\parallel_2$ and dot product $\langle \psi_1 | \psi_2 \rangle$ of states can easily be computed using the [norm](documentation/algebra/algebra.md/#norm) and [dot/dotC](documentation/algebra/algebra.md/#dot) functions.

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:usage_guide_stat4"
	```
	
=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:usage_guide_stat4"
	```

The function [dotC](documentation/algebra/algebra.md/#dot) is only available in C++, and returns a complex (double precision) number whenever one of the involved states is complex. This is necessary, as the return type of a function must be known at compile time in C++, whereas Julia permits dynamic typing. The coefficients of a given state can be retrieved using the [vector/vectorC](documentation/states/state.md/#vectorvectorc).

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:usage_guide_stat5"
	```
	
=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:usage_guide_stat5"
	```

Again, the function [vectorC](documentation/states/state.md/#vectorvectorc) only exists in the C++ version since the return type needs to be known at compile time. In Julia, the type of the vector is decided at runtime. 

Finally, we can apply an operator [OpSum](documentation/operators/opsum.md) to a state using the [apply](documentation/algebra/apply.md) function.

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:usage_guide_stat6"
	```
	
=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:usage_guide_stat6"
	```

Importantly, if the block of the [State](documentation/states/state.md) object has a well-defined quantum number, for example, a conserved particle number, XDiag will automatically detect the quantum number of the resulting state or report an error if the operator does not have a well-defined quantum number. This could be the case, for example, when applying a raising or lowering operator on a particle number conserving state. The [apply](documentation/algebra/apply.md) function acts on a state without creating a matrix representation of the operator, sometimes referred to as *on-the-fly* matrix application.
