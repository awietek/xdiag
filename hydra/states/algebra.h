#pragma once

#include <lila/all.h>

#include <hydra/common.h>
#include <hydra/utils/logger.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

#include <hydra/blocks/blocks.h>
#include <hydra/blocks/target_block.h>

#include <hydra/parallel/mpi/dot_mpi.h>

#include <hydra/states/state.h>
#include <hydra/states/zero_state.h>
#include <hydra/states/random_state.h>


namespace hydra {

template <class coeff_t, class Block>
coeff_t Dot(State<coeff_t, Block> const &v, State<coeff_t, Block> const &w) {
  if (v.block() != w.block()) {
    Log.err("Error: trying to perform Dot product of distinct blocks");
  }
  if constexpr (mpi::is_mpi_block<Block>) {
    return DotMPI(v.vector(), w.vector());
  } else {
    return Dot(v.vector(), w.vector());
  }
}

template <class coeff_t, class Block>
coeff_t Norm(State<coeff_t, Block> const &v) {
  return std::sqrt(Dot(v, v));
}

template <class coeff_t, class Block>
coeff_t Inner(BondList const &bonds, Couplings const &couplings,
              State<coeff_t, Block> const &v) {
  auto Hv = ZeroState<coeff_t, Block>(v.block());
  Apply(bonds, couplings, v, Hv);
  return Dot(v, Hv);
}

template <class coeff_t, class Block>
coeff_t Inner(Bond const &bond, State<coeff_t, Block> const &v) {
  BondList bonds;
  bonds << bond;
  auto cpl = bond.coupling();
  Couplings cpls;
  cpls[cpl] = 1.0;
  return Inner(bonds, cpls, v);
}

template <class coeff_t, class Block>
coeff_t Inner(Bond const &bond, complex coeff, State<coeff_t, Block> const &v) {
  BondList bonds;
  bonds << bond;
  auto cpl = bond.coupling();
  Couplings cpls;
  cpls[cpl] = coeff;
  return Inner(bonds, cpls, v);
}

template <class coeff_t, class Block>
void Apply(BondList const &bonds, Couplings const &couplings,
           State<coeff_t, Block> const &state_in,
           State<coeff_t, Block> &state_out) {
  Apply(bonds, couplings, state_in.block(), state_in.vector(),
        state_out.block(), state_out.vector());
}

template <class coeff_t, class Block>
void Apply(Bond const &bond, State<coeff_t, Block> const &state_in,
           State<coeff_t, Block> &state_out) {
  BondList bonds;
  bonds << bond;
  auto cpl = bond.coupling();
  Couplings cpls;
  cpls[cpl] = 1.0;
  Apply(bonds, cpls, state_in, state_out);
}

template <class coeff_t, class Block>
void Apply(Bond const &bond, complex coeff,
           State<coeff_t, Block> const &state_in,
           State<coeff_t, Block> &state_out) {
  BondList bonds;
  bonds << bond;
  auto cpl = bond.coupling();
  Couplings cpls;
  cpls[cpl] = coeff;
  Apply(bonds, cpls, state_in, state_out);
}

template <class coeff_t, class Block>
State<coeff_t, Block> Apply(BondList const &bonds, Couplings const &couplings,
                            State<coeff_t, Block> const &state_in) {
  auto block_out = TargetBlock(bonds, couplings, state_in.block());
  auto state_out = ZeroState<coeff_t>(block_out);
  Apply(bonds, couplings, state_in, state_out);
  return state_out;
}

template <class coeff_t, class Block>
State<coeff_t, Block> Apply(Bond const &bond,
                            State<coeff_t, Block> const &state_in) {
  BondList bonds;
  bonds << bond;
  auto cpl = bond.coupling();
  Couplings cpls;
  cpls[cpl] = 1.0;
  return Apply(bonds, cpls, state_in);
}

template <class coeff_t, class Block>
State<coeff_t, Block> Apply(Bond const &bond, complex coeff,
                            State<coeff_t, Block> const &state_in) {
  BondList bonds;
  bonds << bond;
  auto cpl = bond.coupling();
  Couplings cpls;
  cpls[cpl] = coeff;
  return Apply(bonds, cpls, state_in);
}

} // namespace hydra
