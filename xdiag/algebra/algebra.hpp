#pragma once

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/state.hpp>

namespace xdiag {

// comparison
XDIAG_API bool isapprox(State const &v, State const &w, double rtol = 1e-12,
                        double atol = 1e-12);

// Various norms
XDIAG_API double norm(State const &v);
XDIAG_API double norm1(State const &v);
XDIAG_API double norminf(State const &v);

// dot
XDIAG_API double dot(State const &v, State const &w);
XDIAG_API complex dotC(State const &v, State const &w);

XDIAG_API double inner(OpSum const &ops, State const &v);
XDIAG_API double inner(Op const &op, State const &v);
XDIAG_API complex innerC(OpSum const &ops, State const &v);
XDIAG_API complex innerC(Op const &op, State const &v);

// Internal routines
double dot(Block const &block, arma::vec const &v, arma::vec const &w);
complex dot(Block const &block, arma::cx_vec const &v, arma::cx_vec const &w);

template <typename coeff_t>
double norm(Block const &block, arma::Col<coeff_t> const &v);

template <typename coeff_t>
double norm1(Block const &block, arma::Col<coeff_t> const &v);

template <typename coeff_t>
double norminf(Block const &block, arma::Col<coeff_t> const &v);

// arithmetic operators
XDIAG_API State &operator*=(State &X, double alpha);
XDIAG_API State operator*(State const &X, double alpha);
XDIAG_API State operator*(double alpha, State const &X);

XDIAG_API State &operator*=(State &X, complex alpha);
XDIAG_API State operator*(State const &X, complex alpha);
XDIAG_API State operator*(complex alpha, State const &X);

XDIAG_API State &operator/=(State &X, double alpha);
XDIAG_API State operator/(State const &X, double alpha);

XDIAG_API State &operator/=(State &X, complex alpha);
XDIAG_API State operator/(State const &X, complex alpha);

XDIAG_API State &operator+=(State &v, State const &w);
XDIAG_API State &operator-=(State &v, State const &w);
XDIAG_API State operator+(State const &v, State const &w);
XDIAG_API State operator-(State const &v, State const &w);

} // namespace xdiag
