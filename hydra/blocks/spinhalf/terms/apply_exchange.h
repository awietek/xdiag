#pragma once

#include <functional>

#include <hydra/bitops/bitops.h>
#include <hydra/common.h>

#include <hydra/blocks/spinhalf/terms/apply_term_offdiag_no_sym.h>
#include <hydra/blocks/spinhalf/terms/apply_term_offdiag_sym.h>

namespace hydra::spinhalf {

// Exchange term: J/2 (S^+_i S^-_j + S^-_i S^+_j)

template <typename bit_t, typename coeff_t, bool symmetric, class IndexingIn,
          class IndexingOut, class Fill>
void apply_exchange(Bond const &bond, IndexingIn &&indexing_in,
                    IndexingOut &&indexing_out, Fill &&fill) {

  assert(bond.coupling_defined());
  assert(bond.type_defined() && (bond.type() == "EXCHANGE"));
  assert(bond.size() == 2);
  assert(bond.sites_disjoint());

  coeff_t J = bond.coupling<coeff_t>();
  int64_t s1 = bond[0];
  int64_t s2 = bond[1];
  bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);

  // Define actions of bond
  auto non_zero_term = [&flipmask](bit_t spins) -> bool {
    return bitops::popcnt(spins & flipmask) & 1;
  };
  std::function<std::pair<bit_t, coeff_t>(bit_t)> term_action;

  coeff_t Jhalf = J / 2.0;
  coeff_t Jhalf_conj = hydra::conj(Jhalf);

  if constexpr (is_real<coeff_t>()) {
    (void)Jhalf_conj;
    term_action = [&flipmask,
                   &Jhalf](bit_t spins) -> std::pair<bit_t, coeff_t> {
      bit_t spins_flip = spins ^ flipmask;
      return {spins_flip, Jhalf};
    };
  } else {
    term_action = [&flipmask, &s1, &Jhalf,
                   &Jhalf_conj](bit_t spins) -> std::pair<bit_t, coeff_t> {
      bit_t spins_flip = spins ^ flipmask;
      return {spins_flip, bitops::gbit(spins, s1) ? Jhalf : Jhalf_conj};
    };
  }

  // Dispatch either symmetric of unsymmetric term application
  if constexpr (symmetric) {
    spinhalf::apply_term_offdiag_sym<bit_t, coeff_t>(
        indexing_in, indexing_out, non_zero_term, term_action, fill);
  } else {
    spinhalf::apply_term_offdiag_no_sym<bit_t, coeff_t>(
        indexing_in, indexing_out, non_zero_term, term_action, fill);
  }
}

} // namespace hydra::spinhalf
