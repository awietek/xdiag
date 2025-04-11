# $S=\frac12$ chain: Symmetries and Ground State Correlators

**Author:** Paul Ebert

This example demonstrates how ground state expectation values can be computed for the spin-$\frac{1}{2}$ XXX antiferromagnetic chain

$$
    \mathcal{H} = \sum_{\langle i, j \rangle} \bm{S}_i \cdot \bm{S}_j = \sum_{\langle i, j \rangle} (S^x_i S^x_j + S^y_i S^y_j + S^z_i S^z_j)
$$

with periodic boundary conditions. The Hamiltonian is invariant under translations and decomposes into $N$ irreducible representations (irreps) labelled by lattice momentum $2\pi/N \times k$ where $k=0, 1, \ldots, N-1$. The Lanczos algorithm is run on all momentum sub-blocks, effectively exchanging a longer runtime for less memory consumption. An important caveat when working with representations is that the operators inside expectation values must be symmetrized with respect to the symmetry group.

=== "C++"
	```c++
	--8<-- "examples/spinhalf_chain_gs_corr_symmetries/main.cpp"
	```

=== "Julia"
	```julia
	--8<-- "examples/spinhalf_chain_gs_corr_symmetries/main.jl"
	```
