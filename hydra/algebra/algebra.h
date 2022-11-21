#pragma once

#include <hydra/common.h>
#include <hydra/operators/bondlist.h>
#include <hydra/states/state.h>
#include <hydra/blocks/blocks.h>

namespace hydra {

template <class coeff_t>
coeff_t dot(State<coeff_t> const &v, State<coeff_t> const &w);

template <class coeff_t> real_t<coeff_t> norm(State<coeff_t> const &v);

// inner
template <class coeff_t>
coeff_t inner(BondList const &bonds, State<coeff_t> const &v);
template <class coeff_t>
coeff_t inner(Bond const &bond, State<coeff_t> const &v);
template <class coeff_t>
coeff_t inner(State<coeff_t> const &w, BondList const &bonds,
              State<coeff_t> const &v);
template <class coeff_t>
coeff_t inner(State<coeff_t> const &w, Bond const &bond,
              State<coeff_t> const &v);

// apply
template <class coeff_t>
void apply(BondList const &bonds, Block const &block_in,
           arma::Col<coeff_t> const &vec_in, Block const &block_out,
           arma::Col<coeff_t> &vec_out);
template <class coeff_t>
void apply(BondList const &bonds, State<coeff_t> const &state_in,
           State<coeff_t> &state_out);
template <class coeff_t>
void apply(Bond const &bond, State<coeff_t> const &state_in,
           State<coeff_t> &state_out);

// template <class coeff_t>
// void apply(Bond const &bond, complex coeff, State<coeff_t> const &state_in,
//            State<coeff_t> &state_out);
// template <class coeff_t>
// State<coeff_t> apply(Bond const &bond, State<coeff_t> const &state_in);

// scal
State<complex> &operator/=(State<complex> &X, complex alpha);
State<complex> &operator/=(State<complex> &X, double alpha);
State<double> &operator/=(State<double> &X, double alpha);

State<complex> &operator*=(State<complex> &X, complex alpha);
State<complex> &operator*=(State<complex> &X, double alpha);
State<double> &operator*=(State<double> &X, double alpha);

} // namespace hydra
