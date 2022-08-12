#pragma once

#include <string>

#include <hydra/bitops/bitops.h>
#include <hydra/common.h>

#include <hydra/blocks/spinhalf/terms/apply_term_offdiag_no_sym.h>
#include <hydra/blocks/spinhalf/terms/apply_term_offdiag_sym.h>

namespace hydra::terms::spinhalf {

// S+ or S- term: J S^+_i   OR   J S^-_i

template <typename bit_t, typename coeff_t, bool symmetric, class IndexingIn,
          class IndexingOut, class Fill>
void apply_spsm(Bond const &bond, Couplings const &couplings,
                IndexingIn &&indexing_in, IndexingOut &&indexing_out,
                Fill &&fill, std::string spsm) {
  if (!(bond.size() == 1)) {
    Log.err("Error in spinhalf::apply_ising: bond has {} sites, expected 1",
            bond.size());
  }
  int s = bond[0];
  bit_t mask = ((bit_t)1 << s);
  std::string cpl = bond.coupling();
  coeff_t J = utils::get_coupling<coeff_t>(couplings, cpl);

  if (!lila::close(J, 0.)) {

    if (spsm == "S+") {
      auto non_zero_term = [&mask](bit_t spins) -> bool {
        return !(spins & mask);
      };
      auto term_action = [&mask, &J](bit_t spins) -> std::pair<bit_t, coeff_t> {
        bit_t spins_flip = spins | mask;
        return {spins_flip, J};
      };

      if constexpr (symmetric) {
        apply_term_offdiag_sym<bit_t, coeff_t>(
            indexing_in, indexing_out, non_zero_term, term_action, fill);
      } else {
        apply_term_offdiag_no_sym<bit_t, coeff_t>(
            indexing_in, indexing_out, non_zero_term, term_action, fill);
      }
    } else if (spsm == "S-") {

      auto non_zero_term = [&mask](bit_t spins) -> bool {
        return spins & mask;
      };
      auto term_action = [&mask, &J](bit_t spins) -> std::pair<bit_t, coeff_t> {
        bit_t spins_flip = spins ^ mask;
        return {spins_flip, J};
      };

      if constexpr (symmetric) {
        apply_term_offdiag_sym<bit_t, coeff_t>(
            indexing_in, indexing_out, non_zero_term, term_action, fill);
      } else {
        apply_term_offdiag_no_sym<bit_t, coeff_t>(
            indexing_in, indexing_out, non_zero_term, term_action, fill);
      }

    } else {
      Log.err("Error in spinhalf::apply_spsm: spsm must either be S+ or S-");
    }
  }
}

} // namespace hydra::terms::spinhalf
