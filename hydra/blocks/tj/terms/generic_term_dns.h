#pragma once

#include <functional>

namespace hydra::tj {

template <typename bit_t, typename coeff_t, bool symmetric, bool fermi_ups,
          class IndexingIn, class IndexingOut, class NonZeroTermUps,
          class NonZeroTermDns, class TermAction, class Fill>
void generic_term_dns(IndexingIn &&indexing_in, IndexingOut &&indexing_out,
                      NonZeroTermUps &&non_zero_term_ups,
                      NonZeroTermDns &&non_zero_term_dns,
                      TermAction &&term_action, Fill &&fill) {
  int n_sites = indexing_in.n_sites();
  assert(n_sites == indexing_out.n_sites());
  bit_t sitesmask = ((bit_t)1 << n_sites) - 1;

  if constexpr (symmetric) {

    auto const &irrep = indexing_out.irrep();
    std::vector<coeff_t> bloch_factors;
    if constexpr (is_complex<coeff_t>()) {
      bloch_factors = irrep.characters();
    } else {
      bloch_factors = irrep.characters_real();
    }

    for (idx_t idx_up_in = 0; idx_up_in < indexing_in.n_rep_ups();
         ++idx_up_in) {
      bit_t up_in = indexing_in.rep_ups(idx_up_in);

      if (non_zero_term_ups(up_in)) {

        bit_t not_up_in = (~up_in) & sitesmask;
        idx_t up_in_offset = indexing_in.ups_offset(idx_up_in);
        idx_t up_out_offset = indexing_out.ups_offset(idx_up_in);

        auto dncs_in = indexing_in.dns_for_ups_rep(up_in);
        auto syms_in = indexing_in.syms_ups(up_in);
        auto norms_in = indexing_in.norms_for_ups_rep(up_in);

        auto dncs_out = indexing_out.dns_for_ups_rep(up_in);
        auto syms_out = indexing_out.syms_ups(up_in);
        auto norms_out = indexing_out.norms_for_ups_rep(up_in);

        // trivial stabilizer of up_in -> dns have to be deposited
        if (syms_in.size() == 1) {
          idx_t idx_in = up_in_offset;
          for (bit_t dnc_in : dncs_in) {
            idx_t dn_in = bitops::deposit(dnc_in, not_up_in);
            if (non_zero_term_dns(dn_in)) {
              auto [dn_flip, coeff] = term_action(dn_in);
              if ((dn_flip & up_in) == 0) { // tJ constraint
                bit_t dnc_flip = bitops::extract(dn_flip, not_up_in);
                idx_t idx_dnc_flip = indexing_out.dnsc_index(dnc_flip);
                idx_t idx_out = up_out_offset + idx_dnc_flip;

		// Log("NO stabilizer");
		// Log("{};{}", BSTR(up_in), BSTR(dn_in));
		// Log("{};{}", BSTR(up_in), BSTR(dn_flip));
		// Log("not_up_in: {}", BSTR(not_up_in));
		// Log("dnc_flip: {}", BSTR(dnc_flip));
		// Log("coeff: {}", coeff);
		// Log("in: {} out: {}", idx_in, idx_out);
		// Log("");
		
                if constexpr (fermi_ups) {
                  bool fermi_up = (bool)(bitops::popcnt(up_in) & 1);
                  fill(idx_out, idx_in, fermi_up ? -coeff : coeff);
                } else {
                  fill(idx_out, idx_in, coeff);
                }
              }
            }
	    ++idx_in;
          }
        }

        // non-trivial stabilizer of ups -> dns don't have to be deposited
        else {
          idx_t idx_in = up_in_offset;
          idx_t idx_dn_in = 0;
          for (bit_t dn_in : dncs_in) {
            if (non_zero_term_dns(dn_in)) {
              auto [dn_flip, coeff] = term_action(dn_in);
              auto [idx_dn_flip, fermi_dn, sym] =
                  indexing_out.index_dns_fermi_sym(dn_flip, syms_out, dncs_out);

              if (idx_dn_flip != invalid_index) {
                idx_t idx_out = up_out_offset + idx_dn_flip;
                bool fermi_up = indexing_out.fermi_bool_ups(sym, up_in);
                if constexpr (fermi_ups) {
                  fermi_up ^= (bool)(bitops::popcnt(up_in) & 1);
                }
                coeff_t val = coeff * bloch_factors[sym] *
                              norms_out[idx_dn_flip] / norms_in[idx_dn_in];

		// Log("YES stabilizer");
		// Log("{};{}", BSTR(up_in), BSTR(dn_in));
		// Log("{};{}", BSTR(up_in), BSTR(dn_flip));
		// Log("not_up_in: {}", BSTR(not_up_in));
		// Log("coeff: {}", val);
		// Log("in: {} out: {}", idx_in, idx_out);
		// Log("");
	
                fill(idx_out, idx_in, (fermi_up ^ fermi_dn) ? -val : val);
              }
            }
            ++idx_dn_in;
            ++idx_in;
          }
        }
      }
    }

  } else { // if not symmetric

    idx_t size_dncs_in = indexing_in.size_dncs();
    idx_t size_dncs_out = indexing_out.size_dncs();
    auto ups_and_idces = indexing_in.states_indices_ups();
    for (auto [up_in, idx_up_in] : ups_and_idces) {

      if (non_zero_term_ups(up_in)) {

        idx_t up_in_offset = idx_up_in * size_dncs_in;
        idx_t up_out_offset = idx_up_in * size_dncs_out;
        bit_t not_up_in = (~up_in) & sitesmask;

        auto dncs_in = indexing_in.states_dncs(up_in);
        idx_t idx_in = up_in_offset;
        for (bit_t dnc_in : dncs_in) {
          bit_t dn_in = bitops::deposit(dnc_in, not_up_in);
          if (non_zero_term_dns(dn_in)) {
            auto [dn_flip, coeff] = term_action(dn_in);
            if ((up_in & dn_flip) == 0) {
              bit_t dnc_flip = bitops::extract(dn_flip, not_up_in);
              idx_t idx_dn_flip = indexing_out.index_dncs(dnc_flip);
              idx_t idx_out = up_out_offset + idx_dn_flip;
	  
              if constexpr (fermi_ups) {
                bool fermi_up = (bool)(bitops::popcnt(up_in) & 1);
                fill(idx_out, idx_in, fermi_up ? -coeff : coeff);
              } else {
                fill(idx_out, idx_in, coeff);
              }
            }
          }
          ++idx_in;
        }
      }
    }
  }
}

} // namespace hydra::tj
