#pragma once

#include <xdiag/blocks/blocks.h>
#include <xdiag/common.h>
#include <xdiag/operators/bondlist.h>
#include <xdiag/states/state.h>

namespace xdiag {

// Various norms
double norm(State const &v);
double norm1(State const &v);
double norminf(State const &v);

// dot
double dot(State const &v, State const &w);
complex dotC(State const &v, State const &w);

double inner(BondList const &bonds, State const &v);
complex innerC(BondList const &bonds, State const &v);
double inner(Bond const &bonds, State const &v);
complex innerC(Bond const &bonds, State const &v);

double inner(State const &v, BondList const &bonds, State const &w);
complex innerC(State const &v, BondList const &bonds, State const &w);
double inner(State const &v, Bond const &bonds, State const &w);
complex innerC(State const &v, Bond const &bonds, State const &w);

State &operator*=(State &X, complex alpha);
State &operator*=(State &X, double alpha);
State &operator/=(State &X, complex alpha);
State &operator/=(State &X, double alpha);

// // inner
// double inner(Bond const &bond, State const &v);
// double inner(State const &w, BondList const &bonds, State const &v);
// double inner(State const &w, Bond const &bond, State const &v);
// complex innerC(BondList const &bonds, State const &v);
// complex innerC(Bond const &bond, State const &v);
// complex innerC(State const &w, BondList const &bonds, State const &v);
// complex innerC(State const &w, Bond const &bond, State const &v);

// State &operator*=(State &X, complex alpha);
// State &operator*=(State &X, double alpha);

double dot(block_variant_t const &block, arma::vec const &v,
           arma::vec const &w);
complex dot(block_variant_t const &block, arma::cx_vec const &v,
            arma::cx_vec const &w);

template <typename coeff_t>
double norm(block_variant_t const &block, arma::Col<coeff_t> const &v);

template <typename coeff_t>
double norm1(block_variant_t const &block, arma::Col<coeff_t> const &v);

template <typename coeff_t>
double norminf(block_variant_t const &block, arma::Col<coeff_t> const &v);

} // namespace xdiag
