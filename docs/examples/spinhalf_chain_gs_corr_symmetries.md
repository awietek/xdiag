# $S=\frac12$ chain: Symmetries and Ground State Correlators

This example demonstrates how ground state expectation values can be computed for the spin-$\frac{1}{2}$ XXX antiferromagnetic chain

$$
    H = \sum_{\langle i, j \rangle} \bm{S}_i \cdot \bm{S}_j = \sum_{\langle i, j \rangle} (S^x_i S^x_j + S^y_i S^y_j + S^z_i S^z_j)
$$

with periodic boundary conditions.

This is done in three different ways, showcasing how the conservation of the overall number of up-spins as well as translation symmetry can be used to lower memory requirements.

* The function ```gs_correlator_simple()``` is the most basic way of obtaining the ground state and computing a correlator. It does not use any of the aforementioned symmetries of the system and therefore needs to store a full Hilbert space vector. Although being fast, this can cause memory issues.

* The function ```gs_correlator_Sz_sym()``` exploits that $H$ conserves the number of up-spins, decomposing the Hamiltonian into $N+1$ separate blocks with a fixed number $n_{up} = 0, 1, \ldots, N$ of spins pointing up. The Lanczos algorithm is then performed on each of the sub-blocks which may take longer but yields a ground state vector whose size is bounded by the size of the largest sub-block, substantially reducing the memory footprint.

* The function ```gs_correlator_translation_sym()``` exploits that all nearest-neighbor bonds are equal and periodic boundary conditions are used. Thus, the Hamiltonian is invariant under translations and decomposes into $N$ irreducible representations (irreps) labelled by lattice momentum $2\pi/N \times k$ where $k=0, 1, \ldots, N-1$. Again, the Lanczos algorithm is run on all momentum sub-blocks, effectively exchanging a longer runtime for less memory consumption. An important caveat when working with representations is that the operators inside expectation values must be symmetrized with respect to the symmetry group.


=== "C++"
	```c++
	--8<-- "examples/spinhalf_chain_gs_corr_symmetries/main.cpp"
	```

=== "Julia"
	```julia
	--8<-- "examples/spinhalf_chain_gs_corr_symmetries/main.jl"
	```