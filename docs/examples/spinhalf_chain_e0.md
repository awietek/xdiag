# Groundstate energy

**Author:** Alexander Wietek

We compute the ground state energy of a spin S = 1/2 Heisenberg
chain,

$$
    H = \sum_{\langle i, j \rangle} \bm{S}_i \cdot \bm{S}_j 
$$

without using translational symmetries. This minimal example shows how
to set up a model and run a simple Lanczos algorithm.

=== "C++"
	```c++
	--8<-- "examples/spinhalf_chain_e0/main.cpp"
	```

=== "Julia"
	```julia
	--8<-- "examples/spinhalf_chain_e0/main.jl"
	```
	
