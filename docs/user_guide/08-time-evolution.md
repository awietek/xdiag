#### Time evolution

Besides diagonalization, XDiag also offers iterative algorithms to perform (imaginary-) time evolutions of the form,

$$ |\phi(t)\rangle = e^{-iHt} |\psi_0\rangle \quad \text{or} \quad |\eta(t)\rangle = e^{-\tau H} |\psi_0\rangle. $$

The time evolution can be performed by two distinct algorithms. The first is the memory-efficient Lanczos algorithm described by [Hochbruck and Lubich (1997)](https://epubs.siam.org/doi/10.1137/S0036142995280572) which runs a Lanczos algorithm twice, to first compute the tridiagonal matrix and then build the time-evolved state. The second is efficient algorithm proposed in [Expokit](https://dl.acm.org/doi/10.1145/285861.285868). While this algorithm is computationally more efficient and highly accurate, it has higher memory requirements.

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:usage_guide_iter3"
	```
	
=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:usage_guide_iter3"
	```

The algorithm employed can be set using the optional `algorithm` argument, which by default is set to `lanczos` using the memory-efficient algorithm of [Hochbruck and Lubich (1997)](https://epubs.siam.org/doi/10.1137/S0036142995280572). Also, more direct control of both algorithms is provided by the functions [evolve_lanczos](documentation/algorithms/evolve_lanczos.md) and [time_evolve_expokit](documentation/algorithms/time_evolve_expokit.md) which allow more control and return further data such as the tridiagonal matrix of the Lanczos algorithm or error estimates, respectively. 
