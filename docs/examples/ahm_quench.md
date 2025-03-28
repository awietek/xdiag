# Time evolution of quench

In this example, we will first find the ground state for the free electron on a 2D lattice
$$
H = -t \sum_{\langle ij,\sigma \rangle} c_{i\sigma}^\dagger c_{j\sigma}
$$
then quench it with a Hubbard interaction, i.e., we will time evolve it with the full Hubbard Hamiltonian
$$
H = -t \sum_{\langle ij,\sigma \rangle} c_{i\sigma}^\dagger c_{j\sigma} + U\sum_i n_{i\uparrow}n_{i\downarrow}.
$$. 
As an observable, we will study the total double occupancy
$$
    \hat{O} = \sum_i \langle n_{i\uparrow}n_{i\downarrow} \rangle
$$

![Image title](../img/quench.png){ align=center }

We see that the double occupancy initially increases sharply, before going down again and entering an oscillatory behaviour.


=== "C++"
	```c++
	--8<-- "examples/ahm_quench/main.cpp"
	```

=== "Julia"
	```julia
	--8<-- "examples/ahm_quench/main.jl"
	```