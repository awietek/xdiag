#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/common.h>

#include <hydra/operators/operator_utils.h>

#include <hydra/blocks/spinhalf/terms/apply_term_offdiag_no_sym.h>
#include <hydra/blocks/spinhalf/terms/apply_term_offdiag_sym.h>

namespace hydra::terms::spinhalf {

// Exchange term: J/2 (S^+_i S^-_j + S^-_i S^+_j)

template <typename bit_t, typename coeff_t, bool symmetric, class IndexingIn,
          class IndexingOut, class Fill>
void apply_exchange(Bond const &bond, Couplings const &couplings,
                    IndexingIn &&indexing_in, IndexingOut &&indexing_out,
                    Fill &&fill) {
  if (!(bond.size() == 2)) {
    Log.err("Error in spinhalf::apply_exchange: bond has {} sites, expected 2",
            bond.size());
  }

  utils::check_sites_disjoint(bond);
  int s1 = bond[0];
  int s2 = bond[1];
  bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);

  // Determine the J/2 coefficient
  std::string cpl = bond.coupling();
  auto [J, J_conj] = utils::get_coupling_and_conj<coeff_t>(couplings, cpl);
  coeff_t Jhalf = J / 2.;
  coeff_t Jhalf_conj = J_conj / 2.;
  if constexpr (is_real<coeff_t>()) { // just to supress "unused" warning
    assert(Jhalf == Jhalf_conj);
  }

  // function determining if spins are flippable
  auto non_zero_term = [&flipmask](bit_t spins) -> bool {
    return bitops::popcnt(spins & flipmask) & 1;
  };

  // Function to flip spin and get coefficient
  auto term_action = [&flipmask, &s1, &Jhalf,
                      &Jhalf_conj](bit_t spins) -> std::pair<bit_t, coeff_t> {
    bit_t spins_flip = spins ^ flipmask;
    if constexpr (is_complex<coeff_t>()) {
      return {spins_flip, bitops::gbit(spins, s1) ? Jhalf : Jhalf_conj};
    } else {
      return {spins_flip, Jhalf};
    }
  };

  // Dispatch either symmetric of unsymmetric term application
  if (!lila::close(J, 0.)) {
    if constexpr (symmetric) {
      terms::spinhalf::apply_term_offdiag_sym<bit_t, coeff_t>(
          indexing_in, indexing_out, non_zero_term, term_action, fill);
    } else {
      terms::spinhalf::apply_term_offdiag_no_sym<bit_t, coeff_t>(
          indexing_in, indexing_out, non_zero_term, term_action, fill);
    }
  }
}

} // namespace hydra::terms::spinhalf
