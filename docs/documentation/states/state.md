---
title: State
---

A generic state describing a quantum wave function $|\psi \rangle$.

**Sources:** [state.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/states/state.hpp) · [state.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/states/state.cpp)

## Constructors

A state can be constructed in three ways:

1. By only specifying the block. In this case the state is initialized with all coefficients zero.

	=== "Julia"
		```julia
		State(block::Block; real::Bool = true, n_cols::Int64 = 1)
		```
	=== "C++"	
		```c++
		State(Block const &block, bool real = true, int64_t n_cols = 1);
		```

2. By handing a vector of coefficients.

	=== "Julia"
		```julia
		State(block::Block, vec::Vector{Float64})
		State(block::Block, vec::Vector{ComplexF64})
		```
	=== "C++"	
		```c++ 
		State(Block const &block, arma::vec const &vector);
		State(Block const &block, arma::cx_vec const &vector);
		```


3. By handing a matrix whose columns describe several states at once.

	=== "Julia"
		```julia
		State(block::Block, mat::Matrix{Float64})
		State(block::Block, mat::Matrix{ComplexF64})
		```
	=== "C++"	
		```c++ 
		State(Block const &block, arma::mat const &matrix);
		State(Block const &block, arma::cx_mat const &matrix);
		```

The number of rows of the vector and matrix needs do agree with the size of the block.

| Parameter | Description                                                                                    |   |
|:----------|:-----------------------------------------------------------------------------------------------|---|
| block     | The block of a Hilbertspace on which the state is defined                                      |   |
| real      | Flag whether or not the state has real coefficients                                            |   |
| n_cols    | Number of columns of the state (default 1)                                                     |   |
| vector    | A vector containing the coefficients of the state. Must be same size as block.                 |   |
| matrix    | A matrix containing the coefficients of the state. Number of rows must be same as block size . |   |


## Methods

#### nsites

Returns the number of sites of the block the state is defined on.

=== "Julia"
	```julia
	nsites(state::State)::Int64
	```
=== "C++"	
	```c++
	int64_t nsites(State const &s) const
	```

#### isapprox

Returns whether two states are approximately equal.

=== "Julia"
	```julia
	isapprox(v::State, w::State, rtol::Float64, atol::Float64)::Bool
	```
=== "C++"	
	```c++
	bool isapprox(State const &v, State const &w, double rtol = 1e-12,
	              double atol = 1e-12);
	```

#### isreal
Returns whether the state is real.

=== "Julia"
	```julia
	isreal(state::State)::Bool
	```
=== "C++"	
	```c++
	int64_t isreal(State const &s) const;
	```

#### real
Returns the real part of the State.

=== "Julia"
	```julia
    real(state::State)::State
	```
=== "C++"	
	```c++
	State real(State const &s) const;
	```


#### imag
Returns the imaginary part of the State.

=== "Julia"
	```julia
    imag(state::State)::State
	```
=== "C++"	
	```c++
	State imag(State const &s) const;
	```
		
#### make_complex! / make_complex
Turns a real State into a complex State. Does nothing if the state is already complex

=== "Julia"
	```julia
    make_complex!(state::State)
	```
=== "C++"	
	```c++
	void make_complex(State &s);
	```

#### dim
Returns the dimension of the block the state is defined on.

=== "Julia"
	```julia
	dim(block::Spinhalf)::Int64
	```
=== "C++"	
	```c++
	int64_t dim(State const &s) const;
	```

#### size
Returns the `size` of the block (also equal to `nrows`) times the number of columns `ncols`. For distributed blocks the local size of a Block is not the same as the dimension `dim`, which is the overall dimension of the block across all processes.

=== "Julia"
	```julia
	size(s::State)::Int64
	```
=== "C++"	
	```c++
	int64_t size(State const &s);
	```

#### nrows
Returns number of rows of the local storage.

=== "Julia"
	```julia
	nrows(s::State)::Int64
	```
=== "C++"	
	```c++
	int64_t nrows(State const &s);
	```

#### ncols
Returns number of columns.

=== "Julia"
	```julia
	ncols(s::State)::Int64
	```
=== "C++"	
	```c++
	int64_t ncols(State const &s);
	```

#### col
Returns a state created from the n-th column of the storage. Whether or not the storage is copied can be specified by setting the flag "copy".

=== "Julia"
	```julia
	col(s::State, n::Int64 = 1; copy::Bool = true)::State
	```
=== "C++"	
	```c++
	State col(State const &s, int64_t n, bool copy = true);
	```
	
#### vector/vectorC
Returns a vector from the n-th column of the storage. In C++ use "vector"/"vectorC" to either get a real or complex vector.

=== "Julia"
	```julia
	vector(state::State; n::Int64 = 1, copy::Bool=true)
	# no vectorC method in julia
	```
=== "C++"	
	```c++
	arma::vec vector(State const &s, int64_t n = 0, bool copy = true);
	arma::cx_vec vectorC(State const &s, int64_t n = 0, bool copy = true);
	```

#### matrix/matrixC
Returns matrix representing the storage. In C++ use "matrix"/"matrixC" to either get a real or complex matrix.

=== "Julia"
	```julia
	matrix(state::State, copy::Bool=true)
	# no matrixC method in julia
	```
=== "C++"	
	```c++
	arma::vec matrix(State const &s, bool copy = true);
	arma::cx_vec matrixC(State const &s, bool copy = true);
	```

## Usage Example

=== "Julia"
	```julia
	--8<-- "examples/usage_examples/main.jl:state"
	```
=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:state"
	```



