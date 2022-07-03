#pragma once

#include <hydra/blocks/blocks.h>
#include <hydra/operators/operator_qns.h>

namespace hydra {

template <typename bit_t>
Spinhalf<bit_t> TargetBlock(BondList const &bonds, Couplings const &cpls,
                            Spinhalf<bit_t> const &block) {
  int n_sites = block.n_sites();
  int nup_old = block.n_up();
  int nup_new = utils::spinhalf_nup(bonds, cpls, block);
  if (nup_old == nup_new) {
    return block;
  } else {
    if (nup_new == undefined_qn) {
      return Spinhalf<bit_t>(n_sites);
    } else {
      return Spinhalf<bit_t>(n_sites, nup_new);
    }
  }
}

template <typename bit_t>
SpinhalfMPI<bit_t> TargetBlock(BondList const &bonds, Couplings const &cpls,
                               SpinhalfMPI<bit_t> const &block) {
  int n_sites = block.n_sites();
  int nup_old = block.n_up();
  int nup_new = utils::spinhalf_nup(bonds, cpls, block);
  if (nup_old == nup_new) {
    return block;
  } else {
    return SpinhalfMPI<bit_t>(n_sites, nup_new);
  }
}

} // namespace hydra
