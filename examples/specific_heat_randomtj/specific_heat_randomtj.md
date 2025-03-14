# Specific Heat Random $\emph{t-J}$ Model

We use full ED to calculate the specific heat in t-J model with $N$s with all to all random interactions and hoppings,
$$
\mathcal{H} = \frac{1}{\sqrt{N}}\sum_{i\neq j=0}^{N-1} t_{i j} P c^\dagger_{i\alpha} c_{j\alpha} + \frac{1}{\sqrt{N}} \sum_{i< j=0}^{N-1} J_{ij} \boldsymbol{S}_i \cdot \boldsymbol{S}_j,
$$
where $P$ is the projection on non-doubly occupied sites, $\boldsymbol{S}_i=\frac{1}{2}c^\dagger_{i\alpha} \boldsymbol{\sigma}_{\alpha \beta}c_{j\beta}$ is the spin operator on site $i$. Both the hoppings $t_{ij}=t^\ast_{ji}$ and the exchange interation $J_{ij}$ are independent random numbers with zero mean and variance $\bar{t}^2,\bar{J}^2$. This type of system exhibits a transition from a spin glass to a disordered Fermi liquid with increasing doping [[1]](#1). Notably, at the critical value of the doping, $p_c \sim1/3$, the model displays features reminiscent of the criticality observed in SYK models.

The following algorithm builds the $\emph{t-J}$ Hamiltonian computing the specific heat using full ED for each realization of the random couplings. The specific heat is then computed by performing the disorder averaging in the post-procesing.



=== "C++"
	```c++
	--8<-- "examples/specific_heat_randomtj/main.cpp"
	```

=== "Julia"
	```julia
	--8<-- "examples/specific_heat_randomtj/main.jl"
	```
	
## references
<a id="1" href="https://doi.org/10.1103/PhysRevLett.108.240401">[1]</a>
H. Shackleton, A. Wietek, A. Georges, S. Sachdev, Quantum Phase Transition at Nonzero Doping in a Random $\emph{t-J}$ Model. Phys. Rev. Lett. 126, 136602 (2021)