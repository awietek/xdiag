#pragma once

#include <hydra/blocks/blocks.h>
#include <hydra/common.h>
#include <hydra/operators/bondlist.h>
#include <hydra/parallel/mpi/dot_mpi.h>
#include <hydra/states/random_state.h>
#include <hydra/states/state.h>
#include <hydra/states/zero_state.h>
#include <hydra/utils/logger.h>

namespace hydra {

template <class coeff_t, class Block>
inline coeff_t dot(State<coeff_t, Block> const &v,
                   State<coeff_t, Block> const &w) {
  if (v.block() != w.block()) {
    Log.err("Error: trying to perform Dot product of distinct blocks");
  }
  if constexpr (mpi::is_mpi_block<Block>) {
    return DotMPI(v.vector(), w.vector());
  } else {
    return arma::cdot(v.vector(), w.vector());
  }
}

template <class coeff_t, class Block>
inline real_t<coeff_t> norm(State<coeff_t, Block> const &v) {
  return std::abs(std::sqrt(dot(v, v)));
}

template <class coeff_t, class Block>
inline coeff_t inner(BondList const &bonds, State<coeff_t, Block> const &v) {
  auto Hv = zero_state<coeff_t, Block>(v.block());
  apply(bonds, v, Hv);
  return dot(v, Hv);
}

template <class coeff_t, class Block>
inline coeff_t inner(Bond const &bond, State<coeff_t, Block> const &v) {
  BondList bonds;
  bonds << bond;
  return inner(bonds, v);
}

template <class coeff_t, class Block>
inline void apply(BondList const &bonds, State<coeff_t, Block> const &state_in,
                  State<coeff_t, Block> &state_out) {
  apply(bonds, state_in.block(), state_in.vector(), state_out.block(),
        state_out.vector());
}

template <class coeff_t, class Block>
inline void apply(Bond const &bond, State<coeff_t, Block> const &state_in,
                  State<coeff_t, Block> &state_out) {
  BondList bonds;
  bonds << bond;
  apply(bonds, state_in, state_out);
}

template <class coeff_t, class Block>
inline void apply(Bond const &bond, complex coeff,
                  State<coeff_t, Block> const &state_in,
                  State<coeff_t, Block> &state_out) {
  BondList bonds;
  bonds << bond;
  apply(bonds, state_in, state_out);
}

template <class coeff_t, class Block>
inline State<coeff_t, Block> apply(Bond const &bond,
                                   State<coeff_t, Block> const &state_in) {
  BondList bonds;
  bonds << bond;
  return apply(bonds, state_in);
}

template <class Block>
inline State<complex, Block> &operator/=(State<complex, Block> &X,
                                         complex alpha) {
  X.vector() /= alpha;
  return X;
}

template <class Block>
inline State<complex, Block> &operator/=(State<complex, Block> &X,
                                         double alpha) {
  X.vector() /= alpha;
  return X;
}

template <class Block>
inline State<double, Block> &operator/=(State<double, Block> &X, double alpha) {
  X.vector() /= alpha;
  return X;
}

} // namespace hydra
