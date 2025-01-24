---
title: State
---

A generic state describing a quantum wave function $|\psi \rangle$.

**Sources** [state.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/states/state.hpp) [state.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/states/state.cpp)

## Constructors

A state can be constructed in two ways:

**1.** By only specifying the block. In this case the state is initialized with all coefficients zero.

=== "C++"	
	```c++
	State(Block const &block, bool real = true, int64_t n_cols = 1);
	```	
=== "Julia"
	```julia
	State(block::Block; real::Bool = true, n_cols::Int64 = 1)
	```
---
**2.** By handing a vector or matrix of coefficients.

=== "C++"	
	```c++ 
	template <typename block_t, typename coeff_t>
	State(block_t const &block, arma::Col<coeff_t> const &vector);

	template <typename block_t, typename coeff_t>
	State(block_t const &block, arma::Mat<coeff_t> const &matrix);
	```
=== "Julia"
	```julia
	State(block::Block, vec::Vector{Float64})
	State(block::Block, vec::Vector{ComplexF64})
	State(block::Block, mat::Matrix{Float64})
	State(block::Block, mat::Matrix{ComplexF64})
	```	

| Parameter | Description                                                                                    |   |
|:----------|:-----------------------------------------------------------------------------------------------|---|
| block     | The block of a Hilbertspace on which the state is defined                                      |   |
| real      | Flag whether or not the state has real coefficients                                            |   |
| n_cols    | Number of columns of the state (default 1)                                                     |   |
| vector    | A vector containing the coefficients of the state. Must be same size as block.                 |   |
| matrix    | A matrix containing the coefficients of the state. Number of rows must be same as block size . |   |

---

## Methods

#### nsites

Returns the number of sites of the block the state is defined on.

=== "C++"	
	```c++
	int64_t nsites(State const &s) const
	```
=== "Julia"
	```julia
	nsites(state::State)
	```
---

#### isreal
Returns whether the state is real.
	
=== "C++"	
	```c++
	int64_t isreal(State const &s) const;
	```

=== "Julia"
	```julia
	isreal(state::State)
	```
---

#### real
Returns whether the real part of the State.

=== "C++"	
	```c++
	State real(State const &s) const;
	```
=== "Julia"
	```julia
    real(state::State)
	```
---

#### imag
Returns whether the imaginary part of the State.

=== "C++"	
	```c++
	State imag(State const &s) const;
	```
=== "Julia"
	```julia
    imag(state::State)
	```

---
		
#### make_complex! / make_complex
Turns a real State into a complex State. Does nothing if the state is already complex

=== "C++"	
	```c++
	void make_complex(State &s);
	```
=== "Julia"
	```julia
    make_complex!(state::State)
	```
---

#### dim
Returns the dimension of the block the state is defined on.

=== "C++"	
	```c++
	int64_t dim(State const &s) const;
	```

=== "Julia"
	```julia
	dim(block::Spinhalf)
	```
---

#### size
Returns the size of the block the state is defined on. locally. Same as "dim" for non-distributed Blocks but different for distributed blocks.

=== "C++"	
	```c++
	int64_t size(State const &s);
	```

=== "Julia"
	```julia
	size(s::State)
	```
---

#### nrows
Returns number of rows of the local storage. Same as "size"

=== "C++"	
	```c++
	int64_t nrows(State const &s);
	```

=== "Julia"
	```julia
	nrows(s::State)
	```

---

#### n_cols
Returns number of columns of the local storage.

=== "C++"	
	```c++
	int64_t ncols(State const &s);
	```

=== "Julia"
	```julia
	ncols(s::State)
	```
---

#### col
Returns a state created from the n-th column of the storage. Whether or not the storage is copied can be specified by setting the flag "copy".

=== "C++"	
	```c++
	State col(State const &s, int64_t n, bool copy = true);
	```
	
=== "Julia"
	```julia
	col(s::State, n::Int64 = 1; copy::Bool = true)
	```
---

#### vector/vectorC
Returns a vector from the n-th column of the storage. In C++ use "vector"/"vectorC" to either get a real or complex vector.

=== "C++"	
	```c++
	arma::vec vector(State const &s, int64_t n = 0, bool copy = true);
	arma::cx_vec vectorC(State const &s, int64_t n = 0, bool copy = true);
	```
=== "Julia"
	```julia
	vector(state::State; n::Int64 = 1, copy::Bool=true)
	# no vectorC method in julia
	```
---

#### matrix/matrixC
Returns matrix representing the storage. In C++ use "matrix"/"matrixC" to either get a real or complex matrix.

=== "C++"	
	```c++
	arma::vec matrix(State const &s, bool copy = true);
	arma::cx_vec matrixC(State const &s, bool copy = true);
	```

=== "Julia"
	```julia
	matrix(state::State, copy::Bool=true)
	# no matrixC method in julia
	```
---

## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:state"
	```

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:state"
	```


