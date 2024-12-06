---
title: Op
---

Object describing a single linear operator acting on a Hilbert space. Every operator is defined by a type as a string. Operators can live on a certain set of sites, which can be specified. Also, sometimes the precise form of an operator can be determined by a matrix, which can be specified.

**Source** [op.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/operators/op.hpp)

## Constructors

=== "C++"	
	```c++
	Op() = default;
	explicit Op(std::string type);
	Op(std::string type, int64_t site);
	Op(std::string type, std::vector<int64_t> const &sites);
	Op(std::string type, arma::mat const &matrix);
	Op(std::string type, int64_t site, arma::mat const &matrix);
	Op(std::string type, std::vector<int64_t> const &sites, arma::mat const &matrix);
	Op(std::string type, arma::cx_mat const &matrix);
	Op(std::string type, int64_t site, arma::cx_mat const &matrix);
	Op(std::string type, std::vector<int64_t> const &sites, arma::cx_mat const &matrix);
	```

| Parameter | Description                                                   |          |
|:----------|:--------------------------------------------------------------|----------|
| type      | a string which denotes what kind of operator is represented   |          |
| sites     | defines on which site(s) of the lattice the operator acts on. | optional |
| matrix    | defines a matrix which may be needed to describe an operator. | optional |

!!! warning "1-indexing in Julia / 0-indexing in C++"

	To enumerate the sites of an Op, we start counting at 1 in Julia and 0 in C++.

