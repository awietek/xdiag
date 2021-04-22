#include "electron.h"

#include <cassert>
#include <tuple>

#include <hydra/combinatorics/binomial.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/combinatorics/subsets.h>

namespace hydra {
namespace detail {

template <class bit_t, class SymmetryGroup>
std::tuple<bit_t, bit_t>
representative_electron(bit_t ups, bit_t dns, SymmetryGroup &&symmetry_group) {
  auto [rep_ups, n_sym, sym_ptr] = symmetry_group.representative_indices(ups);
  bit_t rep_dns = std::numeric_limits<bit_t>::max();
  for (int i = 0; i < n_sym; ++i) {
    bit_t tdns = symmetry_group.apply(sym_ptr[i], dns);
    if (tdns < rep_dns)
      rep_dns = tdns;
  }
  return {rep_ups, rep_dns};
}

template <class bit_t, class BitsGen, class SymmetryGroup>
void fill_states_norms_electron(BitsGen &&up_generator, BitsGen &&dn_generator,
                                SymmetryGroup &&symmetry_group,
                                Representation const &irrep,
                                std::vector<bit_t> &upspins,
                                std::vector<bit_t> &dnspins,
                                std::vector<double> &norms) {
  for (auto ups : up_generator) {
    for (auto dns : dn_generator) {
      auto [rep_ups, rep_dns] =
          representative_electron(ups, dns, symmetry_group);

      if ((rep_ups == ups) && (rep_dns == dns)) {
        // Compute norm
        complex amplitude = 0.0;
        for (int sym = 0; sym < (int)irrep.characters.size(); ++sym) {
          bit_t tups = symmetry_group.apply(sym, ups);
          if (tups == rep_ups) {
            bit_t tdns = symmetry_group.apply(sym, dns);
            double fermi_sign_ups = symmetry_group.fermi_sign(sym, ups);
            if (tdns == rep_dns) {
              double fermi_sign_dns = symmetry_group.fermi_sign(sym, dns);
              amplitude +=
                  fermi_sign_ups * fermi_sign_dns * irrep.characters[sym];
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
}

void check_electron_nup_ndn(int n_sites, int nup, int ndn) {
  if ((nup < 0) || (nup > n_sites))
    HydraLog.err("Error creating Electron: "
                 "invalid value of nup");
  if ((ndn < 0) || (ndn > n_sites))
    HydraLog.err("Error creating Electron: "
                 "invalid value of ndn");
}

} // namespace detail

template <class bit_t, class SymmetryGroup>
Electron<bit_t, SymmetryGroup>::Electron(int n_sites, int nup, int ndn)
    : n_sites_(n_sites), charge_conserved_(true), charge_(nup + ndn),
      sz_conserved_(true), sz_(nup - ndn), nup_(nup), ndn_(ndn),
      lintable_ups_(n_sites, nup_), lintable_dns_(n_sites, ndn_),
      symmetry_group_defined_(false),
      size_(combinatorics::binomial(n_sites, nup_) *
            combinatorics::binomial(n_sites, ndn_)) {
  detail::check_electron_nup_ndn(n_sites, nup, ndn);

  bit_t b_up = ((bit_t)1 << nup_) - 1;
  bit_t b_dn = ((bit_t)1 << ndn_) - 1;
  bit_t e_up = combinatorics::get_next_pattern(b_up << (n_sites - nup_));
  bit_t e_dn = combinatorics::get_next_pattern(b_up << (n_sites - ndn_));

  auto advance = [this](bit_t &ups, bit_t &dns, idx_t &idx) {
    dns = combinatorics::get_next_pattern(dns);
    ++idx;
    if (dns >= ((bit_t)1 << this->n_sites_)) {
      dns = ((bit_t)1 << this->ndn_) - 1;
      ups = combinatorics::get_next_pattern(ups);
    }
  };

  begin_ = ElectronIterator<bit_t>(b_up, b_dn, 0, advance);
  end_ = ElectronIterator<bit_t>(e_up, e_dn, size_, advance);
}

template <class bit_t, class SymmetryGroup>
Electron<bit_t, SymmetryGroup>::Electron(int n_sites, int nup, int ndn,
                                         SymmetryGroup symmetry_group,
                                         Representation irrep)
    : n_sites_(n_sites), charge_conserved_(true), charge_(nup + ndn),
      sz_conserved_(true), sz_(nup - ndn), nup_(nup), ndn_(ndn),
      lintable_ups_(n_sites, nup_), lintable_dns_(n_sites, ndn_),
      symmetry_group_defined_(true), symmetry_group_(symmetry_group),
      irrep_(irrep) {
  detail::check_electron_nup_ndn(n_sites, nup, ndn);
  detail::fill_states_norms_electron(
      Combinations<bit_t>(n_sites, nup), Combinations<bit_t>(n_sites, ndn),
      symmetry_group, irrep, upspins_, dnspins_, norms_);

  size_ = (idx_t)upspins_.size();
  assert((idx_t)dnspins_.size() == size_);
  assert((idx_t)norms_.size() == size_);
  auto advance = [this](bit_t &ups, bit_t &dns, idx_t &idx) {
    ups = this->upspins_[++idx];
    dns = this->dnspins_[idx];
  };
  bit_t b_up = (size_ > 0) ? upspins_[0] : 0;
  bit_t b_dn = (size_ > 0) ? dnspins_[0] : 0;
  bit_t e_up = (size_ > 0) ? upspins_[size_ - 1] : 0;
  bit_t e_dn = (size_ > 0) ? dnspins_[size_ - 1] : 0;
  begin_ = ElectronIterator<bit_t>(b_up, b_dn, 0, advance);
  end_ = ElectronIterator<bit_t>(e_up, e_dn, size_, advance);
}

template class Electron<uint16, SpaceGroup<uint16>>;
template class Electron<uint32, SpaceGroup<uint32>>;
template class Electron<uint64, SpaceGroup<uint64>>;

} // namespace hydra
