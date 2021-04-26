#include "electron.h"

#include <cassert>
#include <tuple>

#include <hydra/combinatorics/binomial.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/combinatorics/subsets.h>

namespace hydra {

template <class bit_t, class SymmetryGroup>
Electron<bit_t, SymmetryGroup>::Electron(int n_sites, int nup, int ndn,
                                         SymmetryGroup symmetry_group,
                                         Representation irrep)
    : n_sites_(n_sites), charge_conserved_(true), charge_(nup + ndn),
      sz_conserved_(true), sz_(nup - ndn), nup_(nup), ndn_(ndn),
      lintable_ups_(n_sites, nup_), lintable_dns_(n_sites, ndn_),
      symmetric_(true), symmetry_group_(symmetry_group), irrep_(irrep) {

  electron::check_nup_ndn(n_sites, nup, ndn);

  for (auto ups : Combinations(n_sites, nup)) {
    for (auto dns : Combinations(n_sites, ndn)) {
      auto [rep_ups, rep_dns] = representative(ups, dns);

      if ((rep_ups == ups) && (rep_dns == dns)) {
        // Compute norm
        complex amplitude = 0.0;
        for (int sym = 0; sym < (int)irrep.size(); ++sym) {
          bit_t tups = symmetry_group.apply(sym, ups);
          if (tups == rep_ups) {
            bit_t tdns = symmetry_group.apply(sym, dns);
            double fermi_sign_ups = symmetry_group.fermi_sign(sym, ups);
            if (tdns == rep_dns) {
              double fermi_sign_dns = symmetry_group.fermi_sign(sym, dns);
              amplitude +=
                  fermi_sign_ups * fermi_sign_dns * irrep.character(sym);
            }
          }
        }

        // Append states
        if (std::abs(amplitude) > 1e-8) {
          upspins.push_back(rep_ups);
          dnspins.push_back(rep_dns);
          norms.push_back(std::sqrt(std::abs(amplitude)));
        }
      }
    }
  }

  size_ = (idx_t)upspins_.size();
}

template <class bit_t, class SymmetryGroup>
std::tuple<bit_t, bit_t>
Electron<bit_t, SymmetryGroup>::representative(bit_t ups, bit_t dns) {
  auto [rep_ups, n_sym, sym_ptr] = symmetry_group_.representative_indices(ups);
  bit_t rep_dns = std::numeric_limits<bit_t>::max();
  for (int i = 0; i < n_sym; ++i) {
    bit_t tdns = symmetry_group.apply(sym_ptr[i], dns);
    if (tdns < rep_dns)
      rep_dns = tdns;
  }
  return {rep_ups, rep_dns};
}

template <class bit_t, class SymmetryGroup>
bool Electron<bit_t, SymmetryGroup>::operator==(
    Electron<bit_t, SymmetryGroup> const &rhs) {
  return (n_sites_ == rhs.n_sites_) &&
         (charge_conserved_ == rhs.charge_conserved_) &&
         (charge_ == rhs.charge_) && (sz_conserved_ == rhs.sz_conserved_) &&
         (sz_ == rhs.sz_) && (nup_ == rhs.nup_) && (ndn_ == rhs.ndn_) &&
         (symmetric_ == rhs.symmetric_) &&
         (symmetry_group_ == rhs.symmetry_group_) && (irrep_ == rhs.irrep_);
}

template <class bit_t, class SymmetryGroup>
bool Electron<bit_t, SymmetryGroup>::operator!=(
    Electron<bit_t, SymmetryGroup> const &rhs) {
  return !operator==(rhs);
}

template class Electron<uint16, SpaceGroup<uint16>>;
template class Electron<uint32, SpaceGroup<uint32>>;
template class Electron<uint64, SpaceGroup<uint64>>;

} // namespace hydra
