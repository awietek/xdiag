#pragma once

#include <hydra/common.h>

namespace hydra::terms::spinhalf {

template <typename bit_t, typename coeff_t, class IndexingIn, class IndexingOut,
          class NonZeroTerm, class TermAction, class Fill>
void apply_term_offdiag_sym(IndexingIn &&indexing_in,
                            IndexingOut &&indexing_out,
                            NonZeroTerm &&non_zero_term,
                            TermAction &&term_action, Fill &&fill) {
  Representation irrep_out = indexing_out.irrep();
  std::vector<coeff_t> characters;
  if constexpr (is_complex<coeff_t>()) {
    characters = irrep_out.characters();
  } else {
    characters = irrep_out.characters_real();
  }

  // Run through all states and apply term
  for (auto [spins_in, idx_in] : indexing_in) {

    if (non_zero_term(spins_in)) {
      auto [spins_out, coeff] = term_action(spins_in);
      auto [idx_out, sym] = indexing_out.index_sym(spins_out);

      if (idx_out != invalid_index) {
        double norm_out = indexing_out.norm(idx_out);
        double norm_in = indexing_in.norm(idx_in);
        coeff_t bloch = characters[sym];
        coeff_t val = coeff * bloch * norm_out / norm_in;
        fill(idx_out, idx_in, val);
      }
    }
  }
}

} // namespace hydra::terms::spinhalf
