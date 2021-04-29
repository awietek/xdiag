#include "electron_symmetric_matrix.h"

#include <hydra/utils/bitops.h>

#include <hydra/models/electron/electron_utils.h>
#include <hydra/models/electron/terms/electron_u.h>

namespace hydra {

template <class bit_t, class SymmetryGroup>
lila::Matrix<complex>
matrix_cplx(BondList const &bonds, Couplings const &couplings,
            ElectronSymmetric<bit_t, SymmetryGroup> const &block_in,
            ElectronSymmetric<bit_t, SymmetryGroup> const &block_out) {
  using utils::gbit;
  using utils::popcnt;

  assert(block_in == block_out); // only temporary
  idx_t dim = block_in.size();

  auto symmetry_group = block_in.symmetry_group();
  auto irrep = block_in.irrep();

  if (block_in.size() == 0)
    return lila::Matrix<complex>();

  auto mat = lila::Zeros<complex>(dim, dim);

  // Hubbard U
  if (couplings.defined("U")) {

    if (!couplings.is_real("U")) {
      HydraLog.err("Error creating Electron matrix: "
                   "Hubbard U must be a real number");
    }
    double U = couplings.real("U");
    if (!lila::close(U, 0.)) {
      for (auto [up, lower_upper] : block_in.ups_lower_upper_) {
        idx_t lower = lower_upper.first;
        idx_t upper = lower_upper.second;
        for (idx_t idx = lower; idx < upper; ++idx) {
          bit_t dn = block_in.dn(idx);
          double val = U * popcnt(up & dn);
          mat(idx, idx) += val;
        }
      }
    }
  }

  auto hoppings = bonds.bonds_of_type("HOP");
  auto hoppings_up = bonds.bonds_of_type("HOPUP");
  auto hoppings_dn = bonds.bonds_of_type("HOPDN");
  for (auto hop : hoppings + hoppings_up + hoppings_dn) {

    if (hop.size() != 2)
      HydraLog.err("Error computing Electron hopping: "
                   "hoppings must have exactly two sites defined");

    std::string cpl = hop.coupling();
    if (couplings.defined(cpl) && !lila::close(couplings[cpl], (complex)0.)) {
      int s1 = hop.site(0);
      int s2 = hop.site(1);
      bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
      int l = std::min(s1, s2);
      int u = std::max(s1, s2);
      bit_t spacemask = (((bit_t)1 << (u - l - 1)) - 1) << (l + 1);

      // Apply hoppings on dnspins
      if ((hop.type() == "HOP") || (hop.type() == "HOPDN")) {

        for (auto [up, lower_upper] : block_in.ups_lower_upper_) {
          idx_t lower = lower_upper.first;
          idx_t upper = lower_upper.second;
          auto begin_up = block_in.dns_.begin() + lower;
          auto end_up = block_in.dns_.begin() + upper;

          std::vector<int> stable_syms =
              electron::stabilizer_symmetries(up, symmetry_group);
          auto stable_group = symmetry_group.subgroup(stable_syms);
          auto stable_irrep = irrep.subgroup(stable_syms);

          // trivial stabilizer of ups -> no representative lookup -> faster
          if (stable_syms.size() == 1) {
            for (idx_t idx = lower; idx < upper; ++idx) {
              bit_t dn = block_in.dn(idx);

              // If hopping is possible ...
              if (popcnt(dn & flipmask) == 1) {
                bit_t dn_flip = dn ^ flipmask;

                // Look for index os dn flip
                auto it = std::lower_bound(begin_up, end_up, dn_flip);

                // if a state has been found
                if ((it != end_up) && (*it == dn_flip)) {
                  idx_t idx_out = std::distance(begin_up, it) + lower;

                  complex t = (gbit(dn, s1)) ? couplings[cpl]
                                             : lila::conj(couplings[cpl]);
                  double fermi_hop = popcnt(dn & spacemask) & 1 ? -1. : 1.;
                  complex val = -t * fermi_hop * block_in.norm(idx_out) /
                                block_in.norm(idx);
                  mat(idx_out, idx) += val;
                }
              }
            }
          } else { // non-trivial stabilizer of ups
            for (idx_t idx = lower; idx < upper; ++idx) {
              bit_t dn = block_in.dn(idx);

              // If hopping is possible ...
              if (popcnt(dn & flipmask) == 1) {
                bit_t dn_flip = dn ^ flipmask;

                // Determine the dn representative from stable symmetries
                auto [dn_flip_rep, dn_flip_rep_sym] =
                    stable_group.representative_index(dn_flip);

                // Look for index os dn flip
                auto it = std::lower_bound(begin_up, end_up, dn_flip_rep);

                // if a state has been found
                if ((it != end_up) && (*it == dn_flip_rep)) {
                  idx_t idx_out = std::distance(begin_up, it) + lower;

                  complex t = (gbit(dn, s1)) ? couplings[cpl]
                                             : lila::conj(couplings[cpl]);
                  double fermi_hop = popcnt(dn & spacemask) & 1 ? -1. : 1.;
                  double fermi_up =
                      stable_group.fermi_sign(dn_flip_rep_sym, up);
                  double fermi_dn =
                      stable_group.fermi_sign(dn_flip_rep_sym, dn_flip);
                  complex val = -t * fermi_hop * fermi_up * fermi_dn *
                                stable_irrep.character(dn_flip_rep_sym) *
                                block_in.norm(idx_out) / block_in.norm(idx);

                  mat(idx_out, idx) += val;
                }
              }
            } // non-trivial stabilizer of ups

          } // for (auto [up, lower_upper]
        }
      } //     if ((hop.type() == "HOP") || (hop.type() == "HOPDN")) {

      // Apply hoppings on upspins
      if ((hop.type() == "HOP") || (hop.type() == "HOPUP")) {

        for (auto [dn, lower_upper] : block_in.dns_lower_upper_) {
          idx_t lower = lower_upper.first;
          idx_t upper = lower_upper.second;
          auto begin_dn = block_in.ups_.begin() + lower;
          auto end_dn = block_in.ups_.begin() + upper;

          std::vector<int> stable_syms =
              electron::stabilizer_symmetries(dn, symmetry_group);
          auto stable_group = symmetry_group.subgroup(stable_syms);
          auto stable_irrep = irrep.subgroup(stable_syms);

          // trivial stabilizer of dns -> no representative lookup -> faster
          if (stable_syms.size() == 1) {
            for (idx_t idx = lower; idx < upper; ++idx) {
              bit_t up = block_in.up(idx);

              // If hopping is possible ...
              if (popcnt(up & flipmask) == 1) {
                bit_t up_flip = up ^ flipmask;
                idx_t idx_in = block_in.index_switch_to_index(idx);
                complex chi_switch_in = block_in.character_switch_[idx];

                // Look for index of up_flip
                auto it = std::lower_bound(begin_dn, end_dn, up_flip);

                // if a state has been found
                if ((it != end_dn) && (*it == up_flip)) {
                  idx_t idx_out_switch = std::distance(begin_dn, it) + lower;
                  idx_t idx_out =
                      block_in.index_switch_to_index(idx_out_switch);

                  complex t = (gbit(up, s1)) ? couplings[cpl]
                                             : lila::conj(couplings[cpl]);
                  double fermi_hop = popcnt(up & spacemask) & 1 ? -1. : 1.;
                  complex chi_switch_out =
                      block_in.character_switch_[idx_out_switch];
                  complex val = -t * fermi_hop * lila::conj(chi_switch_in) *
                                chi_switch_out * block_in.norm(idx_out) /
                                block_in.norm(idx_in);
                  mat(idx_out, idx_in) += val;
                }
              }
            }
          } else { // non-trivial stabilizer of ups
            for (idx_t idx = lower; idx < upper; ++idx) {
              bit_t up = block_in.up(idx);

              // If hopping is possible ...
              if (popcnt(up & flipmask) == 1) {
                bit_t up_flip = up ^ flipmask;
                idx_t idx_in = block_in.index_switch_to_index(idx);
                complex chi_switch_in = block_in.character_switch_[idx];

                // Determine the dn representative from stable symmetries
                auto [up_flip_rep, up_flip_rep_sym] =
                    stable_group.representative_index(up_flip);

                // Look for index os dn flip
                auto it = std::lower_bound(begin_dn, end_dn, up_flip_rep);

                // if a state has been found
                if ((it != end_dn) && (*it == up_flip_rep)) {
                  idx_t idx_out_switch = std::distance(begin_dn, it) + lower;
                  idx_t idx_out =
                      block_in.index_switch_to_index(idx_out_switch);

                  complex t = (gbit(up, s1)) ? couplings[cpl]
                                             : lila::conj(couplings[cpl]);
                  double fermi_hop = popcnt(up & spacemask) & 1 ? -1. : 1.;
                  double fermi_up =
                      stable_group.fermi_sign(up_flip_rep_sym, up_flip);
                  double fermi_dn =
                      stable_group.fermi_sign(up_flip_rep_sym, dn);
                  complex chi_switch_out =
                      block_in.character_switch_[idx_out_switch];
                  complex val = -t * fermi_hop * fermi_up * fermi_dn *
                                stable_irrep.character(up_flip_rep_sym) *
                                lila::conj(chi_switch_in) * chi_switch_out *
                                block_in.norm(idx_out) / block_in.norm(idx_in);

                  mat(idx_out, idx_in) += val;
                }
              }
            }
          } // non-trivial stabilizer of ups
        }
      }
    }
  } // for (auto hop : hoppings + hoppings_up + hoppings_dn) {

  return mat;
}

template lila::Matrix<complex> matrix_cplx<uint16, SpaceGroup<uint16>>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint16, SpaceGroup<uint16>> const &block_in,
    ElectronSymmetric<uint16, SpaceGroup<uint16>> const &block_out);
template lila::Matrix<complex> matrix_cplx<uint32, SpaceGroup<uint32>>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint32, SpaceGroup<uint32>> const &block_in,
    ElectronSymmetric<uint32, SpaceGroup<uint32>> const &block_out);
template lila::Matrix<complex> matrix_cplx<uint64, SpaceGroup<uint64>>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint64, SpaceGroup<uint64>> const &block_in,
    ElectronSymmetric<uint64, SpaceGroup<uint64>> const &block_out);

} // namespace hydra
