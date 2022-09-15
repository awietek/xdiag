#pragma once

#include "extern/armadillo/armadillo"

#include <hydra/common.h>

#include <hydra/blocks/blocks.h>
#include <hydra/blocks/spinhalf/spinhalf.h>

#include <hydra/combinatorics/combinations.h>
#include <hydra/combinatorics/subsets.h>

namespace hydra {

template <class bit_t, class coeff_t, class wavefunction_f>
void Fill(Spinhalf<bit_t> const &block, arma::Col<coeff_t> &vec,
          wavefunction_f wavefunction) {
  using combinatorics::Combinations;
  using combinatorics::Subsets;

  assert(block.size() == vec.size());

  int n_sites = block.n_sites();
  int n_up = block.n_up();
  idx_t idx = 0;

  if (block.sz_conserved()) {
    for (bit_t spins : Combinations<bit_t>(n_sites, n_up)) {
      vec(idx++) = wavefunction(spins);
    }
  } else {
    for (bit_t spins : Subsets<bit_t>(n_sites)) {
      vec(idx++) = wavefunction(spins);
    }
  }
}

} // namespace hydra
