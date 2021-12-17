#pragma once

#include <hydra/common.h>

#include <hydra/blocks/utils/block_utils.h>

#include <hydra/bitops/bitops.h>

#include <hydra/indexing/electron/electron_indexing.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/operators/operator_utils.h>

#include <hydra/combinatorics/combinations.h>

namespace hydra::terms {

template <typename bit_t, typename coeff_t, class Indexing, class Filler>
void electron_hopping(BondList const &bonds, Couplings const &couplings,
                      Indexing &&indexing, Filler &&fill) {
  using bitops::gbit;
  using bitops::popcnt;

  idx_t size_up = indexing.size_ups();
  idx_t size_dn = indexing.size_dns();

  auto clean_bonds =
      utils::clean_bondlist(bonds, couplings, {"HOP", "HOPUP", "HOPDN"}, 2);

  for (auto bond : clean_bonds) {

    std::string type = bond.type();
    std::string cpl = bond.coupling();

    utils::check_sites_disjoint(bond);
    int s1 = bond[0];
    int s2 = bond[1];

    // Prepare bitmasks
    bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
    int l = std::min(s1, s2);
    int u = std::max(s1, s2);
    bit_t fermimask = (((bit_t)1 << (u - l - 1)) - 1) << (l + 1);

    // Apply hoppings on upspins
    if ((type == "HOP") || (type == "HOPUP")) {
      idx_t idx_up = 0;
      for (auto up : indexing.states_ups()) {
        if (popcnt(up & flipmask) == 1) {
          coeff_t t;
          if constexpr (is_complex<coeff_t>())
            t = (gbit(up, s1)) ? couplings[cpl] : lila::conj(couplings[cpl]);
          else
            t = lila::real(couplings[cpl]);

          double fermi = popcnt(up & fermimask) & 1 ? -1. : 1.;
          coeff_t val = -t * fermi;

          bit_t up_flip = up ^ flipmask;
          idx_t idx_up_flip = indexing.index_ups(up_flip);
          for (idx_t idx_dn = 0; idx_dn < size_dn; ++idx_dn) {
            idx_t idx_out = idx_up_flip * size_dn + idx_dn;
            idx_t idx_in = idx_up * size_dn + idx_dn;
            fill(idx_out, idx_in, val);
          }
        }
        ++idx_up;
      }
    }

    // Apply hoppings on dnspins
    if ((type == "HOP") || (type == "HOPDN")) {
      idx_t idx_dn = 0;
      for (auto dn : indexing.states_dns()) {
        if (popcnt(dn & flipmask) == 1) {
          coeff_t t;
          if constexpr (is_complex<coeff_t>())
            t = (gbit(dn, s1)) ? couplings[cpl] : lila::conj(couplings[cpl]);
          else
            t = lila::real(couplings[cpl]);

          double fermi = popcnt(dn & fermimask) & 1 ? -1. : 1.;
          coeff_t val = -t * fermi;

          bit_t dn_flip = dn ^ flipmask;
          idx_t idx_dn_flip = indexing.index_dns(dn_flip);

          // TODO: replace this by lila slicing once available
          for (idx_t idx_up = 0; idx_up < size_up; ++idx_up) {
            idx_t idx_out = idx_up * size_dn + idx_dn_flip;
            idx_t idx_in = idx_up * size_dn + idx_dn;
            fill(idx_out, idx_in, val);
          }
        }
        ++idx_dn;
      }
    }
  } // for (auto bond :...
}

} // namespace hydra::terms
