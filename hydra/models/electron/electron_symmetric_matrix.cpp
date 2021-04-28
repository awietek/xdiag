#include "electron_matrix.h"

#include <hydra/combinatorics/combinations.h>
#include <hydra/utils/bitops.h>

#include <hydra/models/electron/terms/electron_hopping.h>
#include <hydra/models/electron/terms/electron_u.h>

namespace hydra {

template <class bit_t>
lila::Matrix<double> matrix_real(BondList const &bonds,
                                 Couplings const &couplings,
                                 ElectronSymmetric<bit_t> const &block_in,
                                 ElectronSymmetric<bit_t> const &block_out) {
  assert(block_in == block_out); // only temporary
  idx_t dim = block_in.size();

  auto mat = lila::Zeros<double>(dim, dim);

  // Hubbard U
  electron::do_U(couplings, block_in,
                 [&mat](idx_t idx, double val) { mat(idx, idx) += val; });

  auto symmetry_group = block_in.symmetry_group();
  auto irrep = block_in.irrep();

  for (auto hop : hoppings + hoppings_up + hoppings_dn) {

    if (hop.size() != 2)
      HydraLog.err("Error computing Electron hopping: "
                   "hoppings must have exactly two sites defined");

    std::string cpl = hop.coupling();
    if (couplings.defined(cpl) && !lila::close(couplings[cpl], (complex)0.)) {
      double t = lila::real(couplings[cpl]);
      int s1 = hop.site(0);
      int s2 = hop.site(1);
      bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
      int l = std::min(s1, s2);
      int u = std::max(s1, s2);
      bit_t spacemask = (((bit_t)1 << (u - l - 1)) - 1) << (l + 1);

      // Apply hoppings on dnspins
      if ((hop.type() == "HOP") || (hop.type() == "HOPDN")) {

        bit_t up = 0;
        auto up_stabilizer_group = symmetry_group;
        auto up_stabilizer_irrep = irrep;

        idx_t idx = 0;

        // Loop over all states
        while (idx < dim) {
          up = block_in.up(idx);

          // Compute new stabilizer group and irrep if upspins change
          std::vector<int> stable_symmetries;
          std::vector<complex> stable_characters;
          for (int sym = 0; sym < symmetry_group.size(); ++sym) {
            bit_t tup = symmetry_group.apply(sym, up);
            assert(tup >= up);
            if (tup == up) {
              stable_symmetries.push_back(sym);
              stable_symmetries.push_back(irrep.character(idx));
            }
          }
          up_stabilizer_group = symmetry_group.subgroup(stable_symmetries);
          up_stabilizer_irrep = Representation(stable_characters);

          if (up_stabilizer_group.size() == 1) {
            while (block_in.up(idx) == up) {
              bit_t dn = block_in.dn(idx);
              if (popcnt(dn & flipmask) == 1) {
                double fermi = popcnt(dn & spacemask) & 1 ? -1. : 1.;
                double val = -t * fermi;
                bit_t dn_flip = dn ^ flipmask;
		

              }
              ++idx;
            }
          }
        }
      }
    }
  }
  return mat;
}

} // namespace hydra
