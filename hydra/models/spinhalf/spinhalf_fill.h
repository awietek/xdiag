#pragma once

#include <lila/all.h>

#include <hydra/common.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/models/models.h>
#include <hydra/models/spinhalf/spinhalf.h>

namespace hydra {

template <class bit_t, class coeff_t, class wavefunction_f>
void Fill(Spinhalf<bit_t> const &block, lila::Vector<coeff_t> &vec,
          wavefunction_f wavefunction) {
  using combinatorics::Combinations;

  assert(block.size() == vec.size());

  int n_sites = block.n_sites();
  int n_up = block.n_up();
  idx_t idx = 0;
  for (bit_t spins : Combinations<bit_t>(n_sites, n_up)) {
    vec(idx++) = wavefunction(spins);
  }
}

} // namespace hydra
