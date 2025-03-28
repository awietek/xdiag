# Tower of State $J_1 - J_2$ Model in the Triangular Lattice

We perform a tower of states (TOS) analysis [[1]](#1) of the $J_1-J_2$ spin-$\frac{1}{2}$ Heisenberg model on the triangular lattice [[2]](#2). This model consists of spin-$\frac{1}{2}$ sites with nearest and next-nearest neighbor Heisenberg interactions,

$$
\mathcal{H} = J_1 \sum_{\langle i, j \rangle} \boldsymbol{S}_i \cdot \boldsymbol{S}_j + J_2 \sum_{\langle\langle i, j \rangle\rangle} \boldsymbol{S}_i \cdot \boldsymbol{S}_j.
$$

The TOS analysis provides strong evidence for spontaneous symmetry breaking (SSB) in the thermodynamic limit, since the ground state of a finite system is completely symmetric. The $J_1-J_2$ Heisenberg model on the triangular lattice can stabilize different types of order depending on the ratio $J_2/J_1$, which break the continuous SO(3) spin rotation symmetry.

To perform the TOS analysis, we converged the lowest-lying eigenvalues using the Lanczos algorithm in each symmetry sector. We then determined the total spin quantum number, $S_{\text{tot}}$, by inspecting, for each energy level, the number of degenerate eigenstates; thus, $S_{\text{tot}}$ is given by the maximum $S^z$. Finally, we plotted the energy spectra as a function of $S_{\text{tot}}\left(S_{\text{tot}} + 1\right)$. 

![Image title](../img/tow_triangular_lattice.png){ align=center }

In the figure above, we show the energy spectra as a function of \(S_{\text{tot}}\left(S_{\text{tot}} + 1\right)\) for \(J_2=0\). In this case, the ground state exhibits a \(120^\circ\) Néel order. As described in [[2]](#2), group representation theory can be used to predict the quantum numbers in the spontaneous symmetry breaking phases. In this cases, the multiplicities of irreducible representations in the Anderson tower of states can be obtained. This can then be confirmed by the exact diagonalization results the plotting script given above also plots the multiplicites for some sectors. We stress that the results presented are for $N_{\text{spins}} = 18$ and, therefore, not all momenta in the first Brillouin zone can be resolved.

=== "C++"
	```c++
	--8<-- "examples/tos_triangular/main.cpp"
	```

=== "Julia"
	```julia
	--8<-- "examples/tos_triangular/main.jl"
	```
	
Plotting Script
=== "Julia"
	```julia
	--8<-- "examples/tos_triangular/plotting.jl"
	```

## references
<a id="1" href="https://journals.aps.org/pr/abstract/10.1103/PhysRev.86.694">[1]</a>
P. W. Anderson, An Approximate Quantum Theory of the Antiferromagnetic Ground State, Phys. Rev. 86, 694 (1952)

<a id="2" href="https://arxiv.org/abs/1704.08622">[3]</a>
Alexander Wietek, Michael Schuler and Andreas M. Läuchli, Studying Continuous Symmetry Breaking using Energy Level Spectroscopy, arXiv:1704.08622 (2017).