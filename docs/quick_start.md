---
title: Quick start
---

## Julia

To install XDiag, enter the package mode in the Julia REPL using `]` and type

```julia
add XDiag
```

That's it! Now we are ready to perform our first exact diagonalization. We compute the ground state energy of the $S=1/2$ Heisenberg chain on a periodic chain lattice in one dimension. The Hamiltonian is given by

$$ H = J\sum_{\langle i,j \rangle} \mathbf{S}_i \cdot \mathbf{S}_j$$

where $\mathbf{S}_i = (S_i^x, S_i^y, S_i^z)$ are the spin $S=1/2$ operators
and $\langle i,j \rangle$ denotes summation over nearest-meighbor sites
$i$ and $j$.

```julia
--8<-- "examples/spinhalf_chain_e0/main.jl"
```

For more examples, see our [Example Collection](examples.md).


## C++ 

Using XDiag with C++, we first need to compile the XDiag library, see the [Library compilation](documentation/compilation/compilation.md#library-compilation) instructions. Then, our application code looks like this:

```c++
--8<-- "examples/spinhalf_chain_e0/main.cpp"
```

The `try / catch` clause implements an error trace mechanism, which we recommend for every XDiag application. The application is then compiled as well, see [Application compilation](documentation/compilation/compilation.md#application-compilation) instructions, resulting in our final executable.

The C++ versions offers the opportunity to optimize compilation the code for the target architecture, see [C++ optimization](documentation/compilation/compilation.md#optimization)
