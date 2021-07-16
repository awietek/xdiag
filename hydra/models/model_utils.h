#pragma once

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

#include <hydra/symmetries/spacegroup.h>
#include <hydra/symmetries/representation.h>

namespace hydra::utils {

bool coupling_is_zero(Bond const &bond, Couplings const &couplings);

template <class bit_t, class SymmetryGroup>
double compute_norm(bit_t ups, bit_t dns, SymmetryGroup &&symmetry_group,
                    Representation const &irrep) {
  assert(symmetry_group.size() == irrep.size());
  complex amplitude = 0.0;
  for (int sym = 0; sym < (int)symmetry_group.size(); ++sym) {
    bit_t tups = symmetry_group.apply(sym, ups);
    if (tups == ups) {
      bit_t tdns = symmetry_group.apply(sym, dns);
      double fermi_sign_ups = symmetry_group.fermi_sign(sym, ups);
      if (tdns == dns) {
        double fermi_sign_dns = symmetry_group.fermi_sign(sym, dns);
        amplitude += fermi_sign_ups * fermi_sign_dns * irrep.character(sym);
      }
    }
  }
  return std::sqrt(std::abs(amplitude));
}

void check_nup_ndn_tj(int n_sites, int nup, int ndn);
void check_nup_ndn_electron(int n_sites, int nup, int ndn);


} // namespace hydra::utils
