#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/common.h>
#include <hydra/operators/bond.h>

namespace hydra::electron {

template <typename bit_t, typename coeff_t, class Indexing, class Filler>
void electron_do_down_flips(bit_t ups, idx_t idx_ups, bit_t flipmask,
                            bit_t fermimask, int64_t sdn, coeff_t val,
                            Indexing &&indexing, Filler &&fill) {
  idx_t size_dns = indexing.size_dns();

  // Get limits of flipped up
  bit_t ups_flip = ups ^ flipmask;
  idx_t idx_ups_flip = indexing.index_ups(ups_flip);

  idx_t idx_out_offset = idx_ups_flip * size_dns;
  idx_t idx_in = idx_ups * size_dns;

  for (auto dns : indexing.states_dns()) {

    if ((bitops::popcnt(dns & flipmask) == 1) &&
        ((bool)bitops::gbit(dns, sdn))) {
      bit_t dns_flip = dns ^ flipmask;
      idx_t idx_out = idx_out_offset + indexing.index_dns(dns_flip);
      fill(idx_out, idx_in, (bitops::popcnt(dns & fermimask) & 1) ? -val : val);
    }

    ++idx_in;
  }
}

template <typename bit_t, typename coeff_t, bool symmetric, class Indexing,
          class Filler>
void apply_exchange(Bond const &bond, Indexing &&indexing, Filler &&fill) {
  assert(bond.coupling_defined());
  assert(bond.type_defined() && (bond.type() == "EXCHANGE"));
  assert(bond.size() == 2);
  assert(bond.sites_disjoint());

  coeff_t J = bond.coupling<coeff_t>();
  int64_t s1 = bond[0];
  int64_t s2 = bond[1];

  // Prepare bitmasks
  bit_t s1mask = (bit_t)1 << s1;
  bit_t s2mask = (bit_t)1 << s2;
  bit_t flipmask = s1mask | s2mask;
  int64_t l = std::min(s1, s2);
  int64_t u = std::max(s1, s2);
  bit_t fermimask = (((bit_t)1 << (u - l - 1)) - 1) << (l + 1);

  if constexpr (symmetric) {

    auto const &irrep = indexing.irrep();

    // Loop over all up configurations
#ifdef _OPENMP
#pragma omp parallel for schedule(guided)
#endif
    for (idx_t idx_ups = 0; idx_ups < indexing.n_rep_ups(); ++idx_ups) {
      bit_t ups = indexing.rep_ups(idx_ups);

      // Exchange gives zero continue
      if ((bitops::popcnt(ups & flipmask) == 2) ||
          (bitops::popcnt(ups & flipmask) == 0)) {
        continue;
      }

      // Set the correct prefactor
      coeff_t Jhalf;
      if constexpr (is_complex<coeff_t>()) {
        Jhalf = (bitops::gbit(ups, s1)) ? J / 2. : hydra::conj(J) / 2.;
      } else {
        Jhalf = J / 2.;
      }

      // Get the flip masks for up and down spins
      bit_t dns_mask = ((ups & flipmask) == s1mask) ? s2mask : s1mask;

      // Compute index and rep of flipped ups
      bit_t ups_flip = ups ^ flipmask;
      idx_t idx_ups_flip = indexing.index_ups(ups_flip);
      bit_t ups_flip_rep = indexing.rep_ups(idx_ups_flip);

      // Get limits, syms, and dns for ingoing ups
      idx_t up_offset_in = indexing.ups_offset(idx_ups);
      auto dnss_in = indexing.dns_for_ups_rep(ups);
      auto norms_in = indexing.norms_for_ups_rep(ups);

      // Get limits, syms, and dns for outgoing ups
      idx_t up_offset_out = indexing.ups_offset(idx_ups_flip);
      auto syms_up_out = indexing.syms_ups(ups_flip);
      auto dnss_out = indexing.dns_for_ups_rep(ups_flip_rep);
      auto norms_out = indexing.norms_for_ups_rep(ups_flip_rep);

      bool fermi_up = !(bool)(bitops::popcnt(ups & fermimask) & 1);

      // trivial up-stabilizer (likely)
      if (syms_up_out.size() == 1) {
        int64_t sym = syms_up_out.front();
        fermi_up ^= indexing.fermi_bool_ups(sym, ups_flip);

        // Fix the bloch factor
        coeff_t prefac;
        if constexpr (is_complex<coeff_t>()) {
          prefac = -Jhalf * irrep.character(sym);
        } else {
          prefac = -Jhalf * real(irrep.character(sym));
        }

        // Fermi-sign of up spins
        bool fermi_up = (bitops::popcnt(ups & fermimask) & 1);
        fermi_up ^= indexing.fermi_bool_ups(sym, ups_flip);

        idx_t idx_dn = 0;
        for (bit_t dns : dnss_in) {

          // If  dns can be raised
          if ((dns & flipmask) == dns_mask) {
            bit_t dns_flip = dns ^ flipmask;
            auto [idx_dn_flip, fermi_dn] =
                indexing.index_dns_fermi(dns_flip, sym, fermimask);

            coeff_t val = prefac / norms_in[idx_dn];
            idx_t idx_in = up_offset_in + idx_dn;
            idx_t idx_out = up_offset_out + idx_dn_flip;
            fill(idx_out, idx_in, (fermi_up ^ fermi_dn) ? -val : val);
          }
          ++idx_dn;
        }
        // non-trivial up-stabilizer (unlikely)
      } else {

        auto syms = syms_up_out;

        // Fix the bloch/prefactors
        std::vector<coeff_t> prefacs(irrep.size());
        if constexpr (is_complex<coeff_t>()) {
          for (int64_t i = 0; i < (int64_t)irrep.size(); ++i) {
            prefacs[i] = -irrep.character(i) * Jhalf;
          }
        } else {
          for (int64_t i = 0; i < (int64_t)irrep.size(); ++i) {
            prefacs[i] = -real(irrep.character(i)) * Jhalf;
          }
        }
        bool fermi_up_hop = (bitops::popcnt(ups & fermimask) & 1);

        idx_t idx_dn = 0;
        for (bit_t dns : dnss_in) {

          // If  dns can be raised
          if ((dns & flipmask) == dns_mask) {
            bit_t dns_flip = dns ^ flipmask;
            auto [idx_dn_flip, fermi_dn, sym] = indexing.index_dns_fermi_sym(
                dns_flip, syms, dnss_out, fermimask);

            if (idx_dn_flip != invalid_index) {

              bool fermi_up =
                  fermi_up_hop ^ indexing.fermi_bool_ups(sym, ups_flip);
              coeff_t val =
                  prefacs[sym] * norms_out[idx_dn_flip] / norms_in[idx_dn];
              idx_t idx_in = up_offset_in + idx_dn;
              idx_t idx_out = up_offset_out + idx_dn_flip;
              fill(idx_out, idx_in, (fermi_up ^ fermi_dn) ? -val : val);
            }
          }
          ++idx_dn;
        }
      } // trivial stabilizer or not
    }   // loop over ups

  } else { // if not symmetric

#ifdef _OPENMP
#pragma omp parallel
    {
      auto ups_and_idces = indexing.states_indices_ups_thread();
#else
    auto ups_and_idces = indexing.states_indices_ups();
#endif

      for (auto [ups, idx_up] : ups_and_idces) {

        if (bitops::popcnt(ups & flipmask) == 1) {

          // Set the correct prefactor
          coeff_t Jhalf;
          if constexpr (is_complex<coeff_t>()) {
            Jhalf = (bitops::gbit(ups, s1)) ? J / 2.0 : hydra::conj(J) / 2.0;
          } else {
            Jhalf = J / 2.;
          }

          // decide Fermi sign of upspins
          if (bitops::popcnt(ups & fermimask) & 1) {
            electron_do_down_flips(ups, idx_up, flipmask, fermimask,
                                   bitops::gbit(ups, s1) ? s2 : s1, Jhalf,
                                   indexing, fill);
          } else {
            electron_do_down_flips(ups, idx_up, flipmask, fermimask,
                                   bitops::gbit(ups, s1) ? s2 : s1, -Jhalf,
                                   indexing, fill);
          }
        }
      }

#ifdef _OPENMP
    }
#endif
  }
}

} // namespace hydra::electron
