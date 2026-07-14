### Hilbert spaces

The first thing to define before any computation is the Hilbert space our model will be defined on. For this, we create an object of type [Spinhalf](documentation/blocks/spinhalf.md) and hand as an argument the number of physical sites $N$.

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:usage_guide_hs1"
	```
=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:usage_guide_hs1"
	```
	
We would like to know which spin configurations the Hilbert space is made up of. To do so, we can iterate over the Hilbert space and print out the configurations and the total Hilbert space dimension.


=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:usage_guide_hs2"
	```
=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:usage_guide_hs2"
	```


XDiag features a convenient way to write logs in C++, with the [Log](documentation/utilities/logging.md) class. The first argument to `Log()` is a format string. In C++ we use the [fmt library](https://fmt.dev/) to be able to write formatted output. The second argument turns our spins into a string. `spins` is of type [ProductState](documentation/states/product_state.md), whose configuration on each site can be individually addressed by using the `[]` operator. The opposite task of computing the index of a given `ProductState` within the block basis can be addressed using the `index` function. An important difference between C++ and Julia is that indices are counted starting from `0` in C++ and `1` in Julia. Hence, also in the above code snippet, the C++ version will start counting the indices from `0` and Julia from `1`.

The precise output depends on which Hilbert space is chosen. At present, XDiag features three distinct types of Hilbert spaces:

* [Spinhalf](documentation/blocks/spinhalf.md): $S=1/2$ spins; each site is either occupied by an $\uparrow$-spin or a $\downarrow$-spin.
* [Electron](documentation/blocks/electron.md): spin $S=1/2$ electrons; each site is either empty $\emptyset$, occupied by an $\uparrow$-spin or $\downarrow$-spin electron, or is doubly occupied $\updownarrow$.
* [tJ](documentation/blocks/tJ.md): spin $S=1/2$ electrons without double occupancies; each site is either empty $\emptyset$, occupied by an $\uparrow$-spin or $\downarrow$-spin electron.

Frequently, many-body systems feature certain symmetries and conservation laws. Common conservation laws include particle number, spin, or momentum conservation. The Hilbert space can then be subdivided into blocks, which are labeled by the respective conserved quantities. Blocks of a Hilbert space with a given particle number can be easily created by handing further arguments when constructing the Hilbert space specifying the particle numbers. The number of $\uparrow$-spins in a [Spinhalf](documentation/blocks/spinhalf.md) block can be specified via,

=== "C++"
	```c++ 
	--8<-- "examples/user_guide/main.cpp:usage_guide_hs3"
	```
=== "Julia"
	```julia 
	--8<-- "examples/user_guide/main.jl:usage_guide_hs3"
	```	

The result of printing out the configurations of specific blocks is shown in the [figure below](#fig-fix-block). This enumeration is important to interpret coefficients of wave functions. By printing out the basis states of spin configurations, the user can also assess how computational basis states are ordered internally in XDiag.


![Fig. 1: Different blocks with fixed up- or down-spins.](img/usage_guide_hs_nosym.png){#fig-fix-block}
