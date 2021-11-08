#pragma once

#include <lila/utils/logger.h>

#include <hydra/bitops/bitops.h>
#include <hydra/blocks/spinhalf/spinhalf.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/common.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

#include <hydra/blocks/utils/block_utils.h>

namespace hydra::terms::spinhalf {

template <class bit_t, class coeff_t, class Filler>
void do_exchange(BondList const &bonds, Couplings const &couplings,
                 Spinhalf<bit_t> const &block, Filler &&fill) {
  using bitops::gbit;
  using combinatorics::Combinations;

  auto exchange = bonds.bonds_of_type("HEISENBERG") +
                  bonds.bonds_of_type("EXCHANGE") + bonds.bonds_of_type("HB");

  for (auto bond : exchange) {

    if (bond.size() != 2)
      lila::Log.err("Error computing Spinhalf Exchange: "
                    "bond must have exactly two sites defined");
    if (utils::coupling_is_non_zero(bond, couplings)) {
      std::string coupling = bond.coupling();
      std::string type = bond.type();
      int s1 = bond.site(0);
      int s2 = bond.site(1);
      if (s1 == s2)
        lila::Log.err("Error computing Spinhalf Exchange: "
                      "operator acting on twice the same site");
      bit_t mask = ((bit_t)1 << s1) | ((bit_t)1 << s2);

      int n_sites = block.n_sites();
      int n_up = block.n_up();

      coeff_t Jhalf;
      coeff_t Jhalf_conj;

      // Set the correct prefactor
      if constexpr (is_complex<coeff_t>()) {
        Jhalf = couplings[coupling] / 2.;
        Jhalf_conj = lila::conj(Jhalf);

      } else {
        Jhalf = lila::real(couplings[coupling] / 2.);
        Jhalf_conj = Jhalf;
      }

      if (lila::close(lila::imag(Jhalf), 0.)) { // Real exchange

	idx_t idx = 0;
        for (bit_t spins : Combinations<bit_t>(n_sites, n_up)) {

          if (bitops::popcnt(spins & mask) & 1) {
            bit_t new_spins = spins ^ mask;
            idx_t new_idx = block.index(new_spins);
            fill(new_idx, idx, Jhalf);
          }
          ++idx;
        }
      } else { // Complex exchange
        if (!is_complex<coeff_t>())
          lila::Log.err(
              "Cannot compute complex exchange with real block in Spinhalf");

	idx_t idx = 0;
        for (bit_t spins : Combinations<bit_t>(n_sites, n_up)) {
          if (bitops::popcnt(spins & mask) & 1) {
            bit_t new_spins = spins ^ mask;
            idx_t new_idx = block.index(new_spins);
            if (gbit(spins, s1)) {
              fill(new_idx, idx, Jhalf);
            } else {
              fill(new_idx, idx, Jhalf_conj);
            }
	  }
	  ++idx;
        }
      }
    }
  }
}

} // namespace hydra::terms::spinhalf
