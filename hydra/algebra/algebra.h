#pragma once

#include <hydra/blocks/blocks.h>
#include <hydra/common.h>
#include <hydra/operators/bondlist.h>
#include <hydra/parallel/mpi/dot_mpi.h>
#include <hydra/states/state.h>
#include <hydra/utils/logger.h>

#include <hydra/blocks/electron/electron_apply.h>
#include <hydra/blocks/spinhalf/spinhalf_apply.h>
#include <hydra/blocks/tj/tj_apply.h>

namespace hydra {

template <class coeff_t>
inline coeff_t dot(State<coeff_t> const &v, State<coeff_t> const &w) {
  if (v.block() != w.block()) {
    Log.err("Error: trying to perform Dot product of distinct blocks");
  }
  if constexpr (mpi::is_mpi_block<Block>) {
    return DotMPI(v.vector(), w.vector());
  } else {
    return arma::cdot(v.vector(), w.vector());
  }
}

template <class coeff_t> inline real_t<coeff_t> norm(State<coeff_t> const &v) {
  return std::abs(std::sqrt(dot(v, v)));
}

template <class coeff_t>
inline coeff_t inner(BondList const &bonds, State<coeff_t> const &v) {
  auto Hv = State<coeff_t>(v.block());
  apply(bonds, v, Hv);
  return dot(v, Hv);
}

template <class coeff_t>
inline coeff_t inner(Bond const &bond, State<coeff_t> const &v) {
  BondList bonds;
  bonds << bond;
  return inner(bonds, v);
}

template <class coeff_t>
inline coeff_t inner(State<coeff_t> const &w, BondList const &bonds,
                     State<coeff_t> const &v) {
  auto Hv = State<coeff_t>(w.block());
  apply(bonds, v, Hv);
  return dot(w, Hv);
}

template <class coeff_t>
inline coeff_t inner(State<coeff_t> const &w, Bond const &bond,
                     State<coeff_t> const &v) {
  BondList bonds;
  bonds << bond;
  return inner(w, bonds, v);
}

template <class coeff_t>
inline void apply(BondList const &bonds, Block const &block_in,
                  arma::Col<coeff_t> const &vec_in, Block const &block_out,
                  arma::Col<coeff_t> &vec_out) {

  std::visit(
      variant::overloaded{
          [&bonds, &vec_in, &vec_out](Spinhalf const &blk_in,
                                      Spinhalf const &blk_out) {
            apply(bonds, blk_in, vec_in, blk_out, vec_out);
          },
          [&bonds, &vec_in, &vec_out](tJ const &blk_in, tJ const &blk_out) {
            apply(bonds, blk_in, vec_in, blk_out, vec_out);
          },
          [&bonds, &vec_in, &vec_out](Electron const &blk_in,
                                      Electron const &blk_out) {
            apply(bonds, blk_in, vec_in, blk_out, vec_out);
          },
          [&bonds, &vec_in, &vec_out](auto const &blk_in, auto const &blk_out) {
            Log.err("Error in apply: Invalid blocks or combination of blocks");
            apply(bonds, blk_in, vec_in, blk_out, vec_out);
          },
      },
      block_in, block_out);
}

template <class coeff_t>
inline void apply(BondList const &bonds, State<coeff_t> const &state_in,
                  State<coeff_t> &state_out) {
  apply(bonds, state_in.block(), state_in.vector(), state_out.block(),
        state_out.vector());
}

template <class coeff_t>
inline void apply(Bond const &bond, State<coeff_t> const &state_in,
                  State<coeff_t> &state_out) {
  BondList bonds;
  bonds << bond;
  apply(bonds, state_in, state_out);
}

template <class coeff_t>
inline void apply(Bond const &bond, complex coeff,
                  State<coeff_t> const &state_in, State<coeff_t> &state_out) {
  BondList bonds;
  bonds << bond;
  apply(bonds, state_in, state_out);
}

template <class coeff_t>
inline State<coeff_t> apply(Bond const &bond, State<coeff_t> const &state_in) {
  BondList bonds;
  bonds << bond;
  return apply(bonds, state_in);
}

template <class Block>
inline State<complex> &operator/=(State<complex> &X, complex alpha) {
  X.vector() /= alpha;
  return X;
}

template <class Block>
inline State<complex> &operator/=(State<complex> &X, double alpha) {
  X.vector() /= alpha;
  return X;
}

template <class Block>
inline State<double> &operator/=(State<double> &X, double alpha) {
  X.vector() /= alpha;
  return X;
}

template <class Block>
inline State<complex> &operator*=(State<complex> &X, complex alpha) {
  X.vector() *= alpha;
  return X;
}

template <class Block>
inline State<complex> &operator*=(State<complex> &X, double alpha) {
  X.vector() *= alpha;
  return X;
}

template <class Block>
inline State<double> &operator*=(State<double> &X, double alpha) {
  X.vector() *= alpha;
  return X;
}

} // namespace hydra
