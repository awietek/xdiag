---
title: State
---

A generic state describing a quantum wave function

**Source** [state.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/states/state.hpp)

## Constructors

=== "Julia"
	```julia
	State(block::Block; real::Bool = true, n_cols::Int64 = 1)
	State(block::Block, vec::Vector{Float64})
	State(block::Block, vec::Vector{ComplexF64})
	State(block::Block, mat::Matrix{Float64})
	State(block::Block, mat::Matrix{ComplexF64})
	```

=== "C++"	
	```c++
    State(Block const &block, bool real = true, int64_t n_cols = 1);
	
	template <typename block_t, typename coeff_t>
	State(block_t const &block, arma::Col<coeff_t> const &vector);
	
	template <typename block_t, typename coeff_t>
	State(block_t const &block, arma::Mat<coeff_t> const &matrix);
	```

| Parameter | Description                                                                                    |   |
|:----------|:-----------------------------------------------------------------------------------------------|---|
| block     | The block of a Hilbertspace on which the state is defined                                      |   |
| real      | Flag whether or not the state has real coefficients                                            |   |
| n_cols    | Number of columns of the state (default 1)                                                     |   |
| vector    | A vector containing the coefficients of the state. Must be same size as block.                 |   |
| matrix    | A matrix containing the coefficients of the state. Number of rows must be same as block size . |   |


## Methods

!!! method "n_sites"

	Returns the number of sites of the block the state is defined on.

	=== "Julia"
		```julia
		n_sites(state::State)
		```

	=== "C++"	
		```c++
		int64_t n_sites() const
		```

!!! method "isreal"
	Returns whether the state is real.
	
	=== "Julia"
		```julia
	    isreal(state::State)
		```

	=== "C++"	
		```c++
		int64_t isreal() const;
		```

!!! method "real"
	Returns whether the real part of the State.
	
	=== "Julia"
		```julia
	    real(state::State)
		```

	=== "C++"	
		```c++
		State real() const;
		```

!!! method "imag"
	Returns whether the imaginary part of the State.
	
	=== "Julia"
		```julia
	    imag(state::State)
		```

	=== "C++"	
		```c++
		State imag() const;
		```
		
!!! method "make_complex! / make_complex"
	Turns a real State into a complex State. Does nothing if the state is already complex
	
	=== "Julia"
		```julia
	    make_complex!(state::State)
		```

	=== "C++"	
		```c++
		void make_complex();
		```


!!! method "dim"
	Returns the dimension of the block the state is defined on.

	=== "Julia"
		```julia
		dim(block::Spinhalf)
		```

	=== "C++"	
		```c++
		int64_t dim() const;
		```


!!! method "size"
	Returns the size of the block  the state is defined on. locally. Same as "dim" for non-distributed Blocks but different for distributed blocks.

	=== "Julia"
		```julia
		size(block::Spinhalf)
		```

	=== "C++"	
		```c++
		int64_t size() const;
		```

!!! method "n_rows"
	Returns number of rows of the local storage. Same as "size"

	=== "Julia"
		```julia
		n_rows(block::Spinhalf)
		```

	=== "C++"	
		```c++
		int64_t n_rows() const;
		```
		
!!! method "n_cols"
	Returns number of columns of the local storage.

	=== "Julia"
		```julia
		n_cols(block::Spinhalf)
		```

	=== "C++"	
		```c++
		int64_t n_cols() const;
		```

!!! method "col"
	Returns a state created from the n-th column of the storage. Whether or not the storage is copied can be specified by setting the flag "copy".

	=== "Julia"
		```julia
		col(state::State, n::Int64 = 1; copy::Bool = true)
		```

	=== "C++"	
		```c++
		State col(int64_t n, bool copy = true) const;
		```

!!! method "vector/vectorC"
	Returns a vector from the n-th column of the storage. In C++ use "vector"/"vectorC" to either get a real or complex vector.

	=== "Julia"
		```julia
		vector(state::State; n::Int64 = 1)
		# no vectorC method in julia
		```

	=== "C++"	
		```c++
		arma::vec vector(int64_t n = 0, bool copy = true) const;
		arma::cx_vec vectorC(int64_t n = 0, bool copy = true) const;
		```

!!! method "matrix/matrixC"
	Returns matrix representing the storage. In C++ use "matrix"/"matrixC" to either get a real or complex matrix.

	=== "Julia"
		```julia
		matrix(state::State)
		# no matrixC method in julia
		```

	=== "C++"	
		```c++
		arma::vec matrix(bool copy = true) const;
		arma::cx_vec matrixC(bool copy = true) const;
		```


## Usage Example

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:state"
	```

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:state"
	```

