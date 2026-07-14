### Measurements

Measurements in the form of expectation values of wavefunctions,

$$   \langle \mathcal{O}\rangle =  \langle \psi | \mathcal{O} | \psi \rangle,$$

can be evaluated using the [inner](documentation/algebra/algebra.md#inner) function. For example, we compute a static spin correlation $\langle S_{0}^{z} S_{j}^{z}\rangle$ between site $0$ (resp. $1$ in Julia) and $j$.

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:usage_guide_measu1"
	```
	
=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:usage_guide_measu1"
	```

Notice, that in Julia sites start counting from $1$, whereas in C++ sites are counted from $0$. Furthermore, if a complex wave function or operator is involved, the function [innerC](documentation/algebra/algebra.md#inner) in C++ should be called, which returns a complex number. In Julia only [inner](documentation/algebra/algebra.md#inner) is available whose return type is decided at runtime. 
