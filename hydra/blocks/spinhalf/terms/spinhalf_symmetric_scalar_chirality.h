#pragma once

#include <hydra/common.h>

#include <hydra/bitops/bitops.h>

#include <hydra/blocks/spinhalf/spinhalf.h>
#include <hydra/blocks/utils/block_utils.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/operators/operator_utils.h>

namespace hydra::terms {

template <typename bit_t, typename coeff_t, class Indexing, class Filler>
void spinhalf_symmetric_scalar_chirality(BondList const &bonds,
                                         Couplings const &couplings,
                                         Indexing &&indexing, Filler &&fill) {
  using bitops::gbit;

  auto bloch_factors = utils::characters<coeff_t>(indexing.irrep());
  auto clean_bonds =
      utils::clean_bondlist(bonds, couplings, {"SCALARCHIRALITY", "SCH"}, 3);
  for (auto bond : clean_bonds) {

    if constexpr (!is_complex<coeff_t>()) {  // coeff_t == real
      lila::Log.warn(
          "Warning: Ignoring SCALARCHIRALITY bonds due to real coefficients");
      break;
    } else {  // coeff_t == complex

      utils::check_sites_disjoint(bond);
      int s1 = bond[0];
      int s2 = bond[1];
      int s3 = bond[2];

      std::string cpl = bond.coupling();
      complex J = utils::get_coupling<coeff_t>(couplings, cpl);
      complex Jquarter = complex(0, 0.25) * J;
      complex Jquarter_conj = lila::conj(Jquarter);

      bit_t spinmask = ((bit_t)1 << s1) | ((bit_t)1 << s2) | ((bit_t)1 << s3);

      // Run through all states and apply bond
      for (idx_t idx_in = 0; idx_in < indexing.size(); ++idx_in) {
        bit_t spins = indexing.state(idx_in);
        bit_t threespins = spins & spinmask;
        bit_t spins_void = spins & (~spinmask);

        // scalar chirality annihilates 000 and 111
        if ((threespins != 0) && (threespins != spinmask)) {
          bit_t b1 = bitops::gbit(spins, s1);
          bit_t b2 = bitops::gbit(spins, s2);
          bit_t b3 = bitops::gbit(spins, s3);
          bit_t threespins_cyclic = (b1 << s2) | (b2 << s3) | (b3 << s1);
          bit_t threespins_acyclic = (b1 << s3) | (b2 << s1) | (b3 << s2);

          bit_t spins_cyclic = spins_void | threespins_cyclic;
          bit_t spins_acyclic = spins_void | threespins_acyclic;

          auto [idx_cyclic, syms_cyclic] = indexing.index_syms(spins_cyclic);
          auto [idx_acyclic, syms_acyclic] = indexing.index_syms(spins_acyclic);

          double norm_in = indexing.norm(idx_in);

          // if new spins has non-zero norm, compute element and fill
          if (idx_cyclic != invalid_index) {
            double norm_cyclic = indexing.norm(idx_cyclic);
            coeff_t bloch = bloch_factors[syms_cyclic[0]];
            coeff_t val = Jquarter * bloch * norm_cyclic / norm_in;
            fill(idx_cyclic, idx_in, val);
          }

          if (idx_acyclic != invalid_index) {
            double norm_acyclic = indexing.norm(idx_acyclic);
            coeff_t bloch = bloch_factors[syms_acyclic[0]];
            coeff_t val = Jquarter_conj * bloch * norm_acyclic / norm_in;
            fill(idx_acyclic, idx_in, val);
          }
        }
      }
    } // if constexpr (!is_complex<coeff_t>())
  }   // for (auto bond : clean_bonds)
}

} // namespace hydra::terms
