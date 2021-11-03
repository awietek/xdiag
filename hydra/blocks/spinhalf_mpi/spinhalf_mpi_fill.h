#pragma once

#include <lila/all.h>

#include <hydra/common.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/blocks/spinhalf_mpi/spinhalf_mpi.h>
#include <hydra/bitops/bitops.h>

namespace hydra {

template <class bit_t, class coeff_t, class wavefunction_f>
void Fill(SpinhalfMPI<bit_t> const &block, lila::Vector<coeff_t> &vec,
          wavefunction_f wavefunction) {
  using combinatorics::Combinations;
  
  assert(block.size() == vec.size());

  int n_postfix_bits = block.n_postfix_bits_;

  idx_t idx=0;
  for (auto prefix : block.prefixes_) {
    int n_up_prefix = bitops::popcnt(prefix);
    int n_up_postfix = block.n_up() - n_up_prefix;
    for (auto postfix : Combinations<bit_t>(n_postfix_bits, n_up_postfix)) {
      bit_t state = (prefix << n_postfix_bits) | postfix;
      vec(idx++) = wavefunction(state);
    }
  }
}

} // namespace hydra
