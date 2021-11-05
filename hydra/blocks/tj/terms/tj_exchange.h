#pragma once

#include <lila/utils/logger.h>

#include <hydra/bitops/bitops.h>
#include <hydra/blocks/tj/tj.h>
#include <hydra/blocks/tj/tj_utils.h>
#include <hydra/blocks/utils/block_utils.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/common.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra::terms::tj {

template <class bit_t, class Filler>
void do_exchange(BondList const &bonds, Couplings const &couplings,
                 tJ<bit_t> const &block, Filler &&fill) {
  using bitops::bits_to_string;
  using bitops::gbit;
  using bitops::popcnt;
  using combinatorics::Combinations;

  int n_sites = block.n_sites();
  int N = n_sites;
  int nup = block.n_up();
  int ndn = block.n_dn();
  int charge = nup + ndn;
  int n_holes = n_sites - charge;
  idx_t size_spins = block.size_spins_;
  auto exchange = bonds.bonds_of_type("HEISENBERG") +
                  bonds.bonds_of_type("EXCHANGE") + bonds.bonds_of_type("HB");
  auto exchange_tj =
      bonds.bonds_of_type("TJHEISENBERG") + bonds.bonds_of_type("TJHB");

  for (auto bond : exchange + exchange_tj) {

    if (bond.size() != 2)
      lila::Log.err("Error computing tJ exchange: "
                    "bonds must have exactly two sites defined");

    if (utils::coupling_is_non_zero(bond, couplings)) {
      std::string cpl = bond.coupling();
      double J = lila::real(couplings[cpl]);
      int s1 = bond[0];
      int s2 = bond[1];
      bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
      bit_t sitesmask = ((bit_t)1 << n_sites) - 1;
      int l = std::min(s1, s2);
      int u = std::max(s1, s2);
      bit_t fermimask = (((bit_t)1 << (u - l - 1)) - 1) << (l + 1);

      idx_t idx = 0;
      idx_t holes_idx = 0;
      for (auto holes : Combinations<bit_t>(n_sites, n_holes)) {
        bit_t not_holes = (~holes) & sitesmask;
        idx_t holes_offset = holes_idx * size_spins;

        for (auto spins : Combinations<bit_t>(charge, nup)) {
          bit_t ups_dns = bitops::deposit(spins, not_holes);

          if (popcnt(ups_dns & flipmask) & 1) { // spins are flippable
            bit_t new_ups_dns = ups_dns ^ flipmask;
            bit_t new_spins = bitops::extract(new_ups_dns, not_holes);

            // lila::Log("s1: {} s2: {}", s1,  s2);
            // lila::Log("holes  : {}", bits_to_string(holes, N));
            // lila::Log("spins  : {}", bits_to_string(spins, N));
            // lila::Log("upsdns : {}", bits_to_string(ups_dns, N));
            // lila::Log("mask   : {}", bits_to_string(flipmask, N));
            // lila::Log("nupsdns: {}", bits_to_string(new_ups_dns, N));
            // lila::Log("nspins : {}\n", bits_to_string(new_spins, N));

	    double fermi_up = (popcnt(ups_dns & fermimask) & 1) ? 1.0 : -1.0;
	    double fermi_dn = (popcnt((~ups_dns) & fermimask) & 1) ? -1.0 : 1.0;
	    double val = J/2 * fermi_up * fermi_dn;

            idx_t new_idx =
                holes_offset + block.lintable_spins_.index(new_spins);
            fill(new_idx, idx, val);
          }
          ++idx;
        }
        ++holes_idx;
      }
    }
  }
}

} // namespace hydra::terms::tj
