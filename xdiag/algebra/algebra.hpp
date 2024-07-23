#pragma once

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/state.hpp>

namespace xdiag {

// Various norms
double norm(State const &v);
double norm1(State const &v);
double norminf(State const &v);

// dot
double dot(State const &v, State const &w);
complex dotC(State const &v, State const &w);

double inner(OpSum const &ops, State const &v);
double inner(Op const &op, State const &v);

complex innerC(OpSum const &ops, State const &v);
complex innerC(Op const &op, State const &v);

double dot(Block const &block, arma::vec const &v, arma::vec const &w);
complex dot(Block const &block, arma::cx_vec const &v, arma::cx_vec const &w);

template <typename coeff_t>
double norm(Block const &block, arma::Col<coeff_t> const &v);

template <typename coeff_t>
double norm1(Block const &block, arma::Col<coeff_t> const &v);

template <typename coeff_t>
double norminf(Block const &block, arma::Col<coeff_t> const &v);

// arithmetic operators
State &operator*=(State &X, double alpha);
State operator*(State const &X, double alpha);
State operator*(double alpha, State const &X);

State &operator*=(State &X, complex alpha);
State operator*(State const &X, complex alpha);
State operator*(complex alpha, State const &X);

State &operator/=(State &X, double alpha);
State operator/(State const &X, double alpha);
  
State &operator/=(State &X, complex alpha);
State operator/(State const &X, complex alpha);

State &operator+=(State &v, State const &w);
State &operator-=(State &v, State const &w);
State operator+(State const &v, State const &w);
State operator-(State const &v, State const &w);

} // namespace xdiag
