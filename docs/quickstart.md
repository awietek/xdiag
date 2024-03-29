---
title: Quick start
---

## Hello World

Let us set up our first program using the `hydra` library. For this we need to create two files. The first is the actual `C++` code.

```C++
#include <hydra/all.h>

using namespace hydra;

int main() try {
  say_hello();
} catch (std::exception const &e) {
  traceback(e);
}
```

The function `say_hello()` prints out a welcome message, which also contains information which exact hydra version is used. What is maybe a bit unfamiliar is the `try / catch` block. Hydra implements a traceback mechanism for runtime errors, which is activated by this idiom. While not stricly necessary here, it is a good practice to make use of this.

Now that the application program is written, we next need to set up the compilation instructions using [CMake](https://cmake.org/). To do so we create a second file called `CMakeLists.txt` in the same directory.

```cmake
cmake_minimum_required(VERSION 3.19)

project(
  hello_world
  LANGUAGES CXX
)

find_package(hydra REQUIRED HINTS "/path/to/hydra/install")
add_executable(main main.cpp)
target_link_libraries(main PUBLIC hydra::hydra)
```

You should replace `"/path/to/hydra/install"` with the appropriate directory where your hydra library is installed after compilation. This exact `CMakeLists.txt` file can be used to compile any hydra application.

!!! info

    For using the distributed hydra library the last line of the above
    `CMakeLists.txt` should be changed to

    ```cmake
    target_link_libraries(main PUBLIC hydra::hydra_distributed)
    ```

We then compile the application code,

```bash
cmake -S . -B build
cmake --build build
```

and finally run our first `hydra` application.

```bash
./build/main
```

## Computing the ground state energy of a spin chain 

We compute the ground state energy of the $S=1/2$ Heisenberg chain on
a periodic chain lattice in one dimension. The Hamiltonian is given by

$$ H = J\sum_{\langle i,j \rangle} \mathbf{S}_i \cdot \mathbf{S}_j$$

where $\mathbf{S}_i = (S_i^x, S_i^y, S_i^z)$ are the spin $S=1/2$ operators
and $\langle i,j \rangle$ denotes summation over nearest-meighbor sites
$i$ and $j$.

The following code, sets up the Hilbert space, defines the Hamiltonian and finally calls an iterative eigenvalue solver to compute the ground state energy.

```C++
#include <hydra/all.h>

using namespace hydra;

int main() try {

  int n_sites = 16;
  int nup = n_sites / 2;

  // Define the Hilbert space block
  auto block = Spinhalf(n_sites, nup);

  // Define the nearest-neighbor Heisenberg Hamiltonian
  BondList bonds;
  for (int i = 0; i < n_sites; ++i) {
    bonds << Bond("HB", "J", {i, (i + 1) % n_sites});
  }

  // Set the coupling constant "J" to one
  bonds["J"] = 1.0;

  // Compute and print the ground state energy
  double e0 = eigval0(bonds, block);
  Log("Ground state energy: {:.12f}", e0);

} catch (std::exception const &e) {
  traceback(e);
}
```