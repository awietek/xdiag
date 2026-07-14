### Operators

Besides Hilbert spaces, the second key objects in quantum mechanics are operators. In a many-body setting, we consider operators of the form,
$$
	O = \sum_{A\subseteq \mathcal{L}} c_A O_A,
$$
where $O_A$ denotes a local operator acting on sites $A=\{a_1, \ldots, a_{l_A}\}$ and $\mathcal{L}$ denotes the lattice and $c_{A}$ are coupling constants. In the case of the Heisenberg model, we would thus have $\mathcal{O}_{A} = \mathbf{S}_i\cdot\mathbf{S}_j$ and $c_A = J$. In XDiag, the local operators are represented via an [Op](documentation/operators/op.md) object while the sum is represented by an [OpSum](documentation/operators/opsum.md) object. These values of the coupling constants $c_A$ can either be a real or complex number, or a string that later needs to be replaced. The Hamiltonian of a spin $S=1 / 2$ Heisenberg chain is created in the following way:

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:usage_guide_op1"
	```
	
=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:usage_guide_op1"
	```

We first create an empty [OpSum](documentation/operators/opsum.md) and then add additional terms to it. The first part of the product denotes the coupling constant, here given as a string. Alternatively, one could have directly used real / complex numbers here. The second part of the product is a single [Op](documentation/operators/op.md) object. It is created with two inputs:

* The type, here `SdotS` denoting an operator of the form $\mathbf{S}_{i} \cdot \mathbf{S}_{j}$. XDiag features a wide variety of different [operator types](documentation/operators/operator_types.md).
* An array defining which site the operator lives on. Notice, that in Julia we start counting the sites from 1, while in C++ we start counting the sites from 0.
