# Tower Of States attractive Hubbard model

We perform a tower of states (TOS) analysis [[1]](#1),[[2]](#2) on the attractive Hubbard model on a square lattice. The (particle-hole symmetric) model can be written as
$$
    \mathcal{H} = -t\sum_{\langle i j \rangle} c_i^{\dagger} c_j - U\sum_{i} (n_{i\uparrow}-\frac{1}{2})(n_{i\downarrow}-\frac{1}{2}).
$$

The attractive Hubbard model is known to display a superconducting ground state, which in mean-field theory is well described by a fermionic U$(1)$ coherent state
$$
    \ket{\psi} = \prod_{k} (u_k + v_k c_{k\uparrow}^\dagger c_{-k\downarrow}^\dagger)\ket{0},
$$
i.e, the ground state is degenerate in the number of particle pairs. We thus expect that in the superfluid state, in the neighbourhood of the most favourable particle number (determined by the chemical potential), there will be a tower of almost degenerate energy states, containing different number of superconducting pairs, but otherwise the same.

To perform the TOS analysis, we gather the eigenvalues using either full exact diagonalization (ED) or the Lanczos algorithm in each symmetry sector, defined by the particle number. Here we use the Lanczos algorithm with $10$ eigenvalues well converged.

![Image title](../img/tos_ahm_Lx(4)_Ly(4).png){ align=center }

As we can see in the spectra above (computed for a $4\times 4$ lattice), for large enough values of $U$, the ground state energy is almost degenerate in the different particle number sectors, providing strong evidence for U$(1)$ spontaneous symmetry breaking.


=== "C++"
	```c++
	--8<-- "examples/tos_ahm/main.cpp"
	```

=== "Julia"
	```julia
	--8<-- "examples/tos_ahm/main.jl"
	```
	
## references
<a id="1" href="https://journals.aps.org/pr/abstract/10.1103/PhysRev.86.694">[1]</a>
P. W. Anderson, An Approximate Quantum Theory of the Antiferromagnetic Ground State, Phys. Rev. 86, 694 (1952)

<a id="2" href="https://arxiv.org/abs/1704.08622">[3]</a>
Alexander Wietek, Michael Schuler and Andreas M. Läuchli, Studying Continuous Symmetry Breaking using Energy Level Spectroscopy, arXiv:1704.08622 (2017).