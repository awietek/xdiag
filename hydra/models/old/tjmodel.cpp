#include <algorithm>
#include <lila/all.h>

#include <hydra/bases/basis_spinhalf.h>
#include <hydra/combinatorics/binomial.h>
#include <hydra/combinatorics/up_down_hole.h>
#include <hydra/utils/bitops.h>

#include "hubbardmodeldetail.h"
#include "tjmodel.h"

namespace hydra {

template <class coeff_t, class bit_t, class idx_t>
TJModel<coeff_t, bit_t, idx_t>::TJModel(BondList bondlist, Couplings couplings,
                                        qn_t qn)
    : n_sites_(bondlist.n_sites()), qn_(qn) {
  using hydra::combinatorics::binomial;

  int nup = qn.n_up;
  int ndn = qn.n_dn;
  assert(nup + ndn <= n_sites_);
  dim_ = binomial(n_sites_, nup) * binomial(n_sites_ - nup, ndn);

  // Currently unused operatuors
  std::vector<std::pair<int, int>> currents_;
  std::vector<coeff_t> current_amplitudes_;
  std::vector<std::pair<int, int>> interactions_;
  std::vector<double> interaction_strengths_;
  double U_;

  detail::set_hubbard_terms<coeff_t>(
      bondlist, couplings, hoppings_, hopping_amplitudes_, currents_,
      current_amplitudes_, interactions_, interaction_strengths_, onsites_,
      onsite_potentials_, szszs_, szsz_amplitudes_, exchanges_,
      exchange_amplitudes_, U_);
}

template <class coeff_t, class bit_t, class idx_t>
lila::Matrix<coeff_t>
TJModel<coeff_t, bit_t, idx_t>::matrix(bool ninj_term) const {

  using hydra::combinatorics::up_hole_to_down;
  using namespace hydra::utils;

  int nup = qn_.n_up;
  int ndn = qn_.n_dn;

  // Try allocating upspin / downspin vectors
  std::vector<bit_t> upspins;
  std::vector<bit_t> dnspins;
  try {
    upspins.resize(dim_);
    dnspins.resize(dim_);
  } catch (...) {
    std::cerr << "Error: Could not allocate upspin/downspin "
              << "vectors for TJModel!" << std::endl
              << std::flush;
    exit(EXIT_FAILURE);
  }

  // Fill vectors holding upspin/downspin configurations
  int64 idx = 0;
  auto basis_up = BasisSpinHalf<bit_t>(n_sites_, nup);
  auto basis_holes_in_up = BasisSpinHalf<bit_t>(n_sites_ - nup, ndn);
  for (auto ups : basis_up)
    for (auto holes : basis_holes_in_up) {
      auto dns = up_hole_to_down(ups.spins, holes.spins);
      upspins[idx] = ups.spins;
      dnspins[idx] = dns;
      ++idx;
    }
  assert(idx == dim_);

  // Define lambda function to get indices
  auto index_of_up_dn = [&upspins, &dnspins](bit_t const &ups,
                                             bit_t const &dns) {
    // Binary search the new indices
    auto up_bounds = std::equal_range(upspins.begin(), upspins.end(), ups);
    auto up_begin = std::distance(upspins.begin(), up_bounds.first);
    auto up_end = std::distance(upspins.begin(), up_bounds.second);
    auto it = std::lower_bound(dnspins.begin() + up_begin,
                               dnspins.begin() + up_end, dns);
    return std::distance(dnspins.begin(), it);
  };

  // Try allocating the matrix
  lila::Matrix<coeff_t> H;
  try {
    H.resize(dim_, dim_);
    Zeros(H);
  } catch (...) {
    std::cerr << "Error: Could not allocate matrix for TJModel!" << std::endl
              << std::flush;
    exit(EXIT_FAILURE);
  }

  // SzSz terms
  int szsz_idx = 0;
  for (auto pair : szszs_) {
    int s1 = pair.first;
    int s2 = pair.second;
    double jz = szsz_amplitudes_[szsz_idx] * 0.25;

    if (std::abs(jz) > 1e-14) {
      // printf("sites %d %d %d %f\n", s1, s2, szsz_idx, jz);

      for (int64 idx = 0; idx < dim_; ++idx) {
        bit_t ups = upspins[idx];
        bit_t dns = dnspins[idx];

        bool up1 = gbit(ups, s1);
        bool up2 = gbit(ups, s2);
        bool dn1 = gbit(dns, s1);
        bool dn2 = gbit(dns, s2);

        if (ninj_term) {
          if ((up1 && up2) || (dn1 && dn2))
            H(idx, idx) += 0;
          else if ((up1 && dn2) || (dn1 && up2))
            H(idx, idx) += -2 * jz;
        } else {
          if ((up1 && up2) || (dn1 && dn2))
            H(idx, idx) += jz;
          else if ((up1 && dn2) || (dn1 && up2))
            H(idx, idx) += -jz;
        }
      }

    } // if (std::abs(jz) > 1e-14)
    ++szsz_idx;
  } // for (auto pair : szszs_)

  // Exchange terms
  int exchange_idx = 0;
  for (auto pair : exchanges_) {
    int s1 = std::min(pair.first, pair.second);
    int s2 = std::max(pair.first, pair.second);
    coeff_t jx = exchange_amplitudes_[exchange_idx] * 0.5;
    bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
    if (std::abs(jx) > 1e-14) {
      for (int64 idx = 0; idx < dim_; ++idx) {
        bit_t ups = upspins[idx];
        bit_t dns = dnspins[idx];
        bool up1 = gbit(ups, s1);
        bool up2 = gbit(ups, s2);
        bool dn1 = gbit(dns, s1);
        bool dn2 = gbit(dns, s2);

        if ((up1 && dn2) || (dn1 && up2)) {
          bit_t flipped_ups = ups ^ flipmask;
          bit_t flipped_dns = dns ^ flipmask;

          // if (popcnt(flipped_ups) != popcnt(ups)) {
          //   using namespace hilbertspaces;
          //   std::cout << PrintBasisSpinHalf(64, ups) << std::endl
          //             << popcnt(ups) << std::endl
          //             << PrintBasisSpinHalf(64, flipped_ups) << std::endl
          //             << popcnt(flipped_ups) << std::endl
          //             << std::endl;
          // }
          // assert(popcnt(flipped_ups) == popcnt(ups));
          // assert(popcnt(flipped_dns) == popcnt(dns));

          int64 flipped_idx = index_of_up_dn(flipped_ups, flipped_dns);
          H(flipped_idx, idx) += jx;
        }

      } // loop over spin configurations
    }   // if (std::abs(jz) > 1e-14)
    ++exchange_idx;
  } // for (auto pair: exchanges_)

  // Hoppings
  int hopping_idx = 0;
  for (auto pair : hoppings_) {
    int s1 = std::min(pair.first, pair.second);
    int s2 = std::max(pair.first, pair.second);
    coeff_t t = hopping_amplitudes_[hopping_idx];
    if (pair.first > pair.second) t = lila::conj(t);
    
    bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
    bit_t firstmask = (bit_t)1 << s1;

    if (std::abs(t) > 1e-14) {
      for (int64 idx = 0; idx < dim_; ++idx) {
        bit_t ups = upspins[idx];
        bit_t dns = dnspins[idx];

        // upspin hopping
        if ((dns & flipmask) == 0) {
          if (((ups & flipmask) != 0) && ((ups & flipmask) != flipmask)) {
            bit_t flipped_ups = ups ^ flipmask;
            int64 flipped_idx = index_of_up_dn(flipped_ups, dns);

            double fermi_up =
                popcnt(gbits(ups ^ dns, s2 - s1, s1)) % 2 == 0 ? 1. : -1.;

            if (ups & firstmask)
              H(flipped_idx, idx) += t * fermi_up;
            else
              H(flipped_idx, idx) -= lila::conj(t) * fermi_up;
          }
        }

        // dnspin hopping
        if ((ups & flipmask) == 0) {
          if (((dns & flipmask) != 0) && ((dns & flipmask) != flipmask)) {
            bit_t flipped_dns = dns ^ flipmask;
            int64 flipped_idx = index_of_up_dn(ups, flipped_dns);

            double fermi_dn =
                popcnt(gbits(ups ^ dns, s2 - s1, s1)) % 2 == 0 ? 1. : -1.;

            if (dns & firstmask)
              H(flipped_idx, idx) += t * fermi_dn;
            else
              H(flipped_idx, idx) -= lila::conj(t) * fermi_dn;
          }
        }
      } // loop over spin configurations
    }   // if (std::abs(t) > 1e-14)
    ++hopping_idx;
  } // for (auto pair : hoppings_)

  return H;
}

template <class coeff_t, class bit_t, class idx_t>
lila::Matrix<coeff_t>
TJModel<coeff_t, bit_t, idx_t>::szMatrix(int siteIndex) const {
  using hydra::combinatorics::up_hole_to_down;
  using namespace hydra::utils;

  assert(0 <= siteIndex &&
         siteIndex < n_sites_); // Check to make sure site index is in range

  int nup = qn_.n_up;
  int ndn = qn_.n_dn;

  // Try allocating upspin / downspin vectors
  std::vector<bit_t> upspins;
  std::vector<bit_t> dnspins;
  try {
    upspins.resize(dim_);
    dnspins.resize(dim_);
  } catch (...) {
    std::cerr << "Error: Could not allocate upspin/downspin "
              << "vectors for TJModel!" << std::endl
              << std::flush;
    exit(EXIT_FAILURE);
  }

  // Fill vectors holding upspin/downspin configurations
  int64 idx = 0;
  auto hs_upspins = BasisSpinHalf<bit_t>(n_sites_, nup);
  auto hs_holes_in_ups = BasisSpinHalf<bit_t>(n_sites_ - nup, ndn);
  for (auto ups : hs_upspins)
    for (auto holes : hs_holes_in_ups) {
      auto dns = up_hole_to_down(ups.spins, holes.spins);
      upspins[idx] = ups.spins;
      dnspins[idx] = dns;
      ++idx;
    }
  assert(idx == dim_);

  // Try allocating the matrix
  lila::Matrix<coeff_t> sz;
  try {
    sz.resize(dim_, dim_);
    Zeros(sz);
  } catch (...) {
    std::cerr << "Error: Could not allocate matrix for TJModel!" << std::endl
              << std::flush;
    exit(EXIT_FAILURE);
  }

  // Iterate through spin configurations

  for (int64 idx = 0; idx < dim_; ++idx) {
    bit_t ups = upspins[idx];
    bit_t dns = dnspins[idx];

    bool up = gbit(ups, siteIndex);
    bool dn = gbit(dns, siteIndex);

    if (up)
      sz(idx, idx) = 0.5;
    else if (dn)
      sz(idx, idx) = -0.5;
    else
      sz(idx, idx) = 0;
  }

  return sz;
}

template <class coeff_t, class bit_t, class idx_t>
lila::Matrix<double>
TJModel<coeff_t, bit_t, idx_t>::sPlusMatrix(int siteIndex) const {
  using hydra::combinatorics::binomial;
  using hydra::combinatorics::up_hole_to_down;
  using namespace hydra::utils;

  assert(0 <= siteIndex &&
         siteIndex < n_sites_); // Check to make sure site index is in range

  int basenup = qn_.n_up;
  int basendn = qn_.n_dn;
  int targetnup = basenup + 1;
  int targetndn = basendn - 1;

  // Make sure target Hilbert space is valid (return a zero matrix otherwise)

  if (targetnup * (targetnup - n_sites_) <= 0 &&
      targetndn * (targetndn - n_sites_) <= 0) {

    int targetdim = binomial(n_sites_, targetnup) *
                    binomial(n_sites_ - targetnup, targetndn);

    // Try allocating upspin / downspin vectors
    std::vector<bit_t> baseupspins;
    std::vector<bit_t> basednspins;
    std::vector<bit_t> targetupspins;
    std::vector<bit_t> targetdnspins;
    try {
      baseupspins.resize(dim_);
      basednspins.resize(dim_);
      targetupspins.resize(targetdim);
      targetdnspins.resize(targetdim);
    } catch (...) {
      std::cerr << "Error: Could not allocate upspin/downspin "
                << "vectors for TJModel!" << std::endl
                << std::flush;
      exit(EXIT_FAILURE);
    }

    // Fill vectors holding upspin/downspin configurations
    int64 idx = 0;
    auto hs_baseupspins = BasisSpinHalf<bit_t>(n_sites_, basenup);
    auto hs_holes_in_baseups =
        BasisSpinHalf<bit_t>(n_sites_ - basenup, basendn);
    for (auto ups : hs_baseupspins)
      for (auto holes : hs_holes_in_baseups) {
        auto dns = up_hole_to_down(ups.spins, holes.spins);
        baseupspins[idx] = ups.spins;
        basednspins[idx] = dns;
        ++idx;
      }
    assert(idx == dim_);

    idx = 0;
    auto hs_targetupspins = BasisSpinHalf<bit_t>(n_sites_, targetnup);
    auto hs_holes_in_targetups =
        BasisSpinHalf<bit_t>(n_sites_ - targetnup, targetndn);
    for (auto ups : hs_targetupspins)
      for (auto holes : hs_holes_in_targetups) {
        auto dns = up_hole_to_down(ups.spins, holes.spins);
        targetupspins[idx] = ups.spins;
        targetdnspins[idx] = dns;
        ++idx;
      }
    assert(idx == targetdim);

    // Try allocating the matrix
    lila::Matrix<double> sPlus;
    try {
      sPlus.resize(targetdim, dim_);
      Zeros(sPlus);
    } catch (...) {
      std::cerr << "Error: Could not allocate matrix for TJModel!" << std::endl
                << std::flush;
      exit(EXIT_FAILURE);
    }

    // Define lambda function to get indices for target space
    auto index_of_up_dn = [&targetupspins, &targetdnspins](bit_t const &ups,
                                                           bit_t const &dns) {
      // Binary search the new indices
      auto up_bounds =
          std::equal_range(targetupspins.begin(), targetupspins.end(), ups);
      auto up_begin = std::distance(targetupspins.begin(), up_bounds.first);
      auto up_end = std::distance(targetupspins.begin(), up_bounds.second);
      auto it = std::lower_bound(targetdnspins.begin() + up_begin,
                                 targetdnspins.begin() + up_end, dns);
      return std::distance(targetdnspins.begin(), it);
    };

    // Iterate through spin configurations
    bit_t flipmask = ((bit_t)1 << siteIndex);

    for (int64 idx = 0; idx < dim_; ++idx) {
      bit_t ups = baseupspins[idx];
      bit_t dns = basednspins[idx];

      bool dn = gbit(dns, siteIndex);

      if (dn) {
        bit_t flipped_ups = ups ^ flipmask;
        bit_t flipped_dns = dns ^ flipmask;
        int64 target_idx = index_of_up_dn(flipped_ups, flipped_dns);
        sPlus(target_idx, idx) = 1;
      }
    }

    return sPlus;

  } else {
    // Return zero matrix

    lila::Matrix<double> sPlus;
    sPlus.resize(1, 1);
    Zeros(sPlus);
    return sPlus;
  }
}

template <class coeff_t, class bit_t, class idx_t>
lila::Matrix<double>
TJModel<coeff_t, bit_t, idx_t>::sMinusMatrix(int siteIndex) const {
  using hydra::combinatorics::binomial;
  using hydra::combinatorics::up_hole_to_down;
  using namespace hydra::utils;

  assert(0 <= siteIndex &&
         siteIndex < n_sites_); // Check to make sure site index is in range

  int basenup = qn_.n_up;
  int basendn = qn_.n_dn;
  int targetnup = basenup - 1;
  int targetndn = basendn + 1;

  if (targetnup * (targetnup - n_sites_) <= 0 &&
      targetndn * (targetndn - n_sites_) <= 0) {

    int targetdim = binomial(n_sites_, targetnup) *
                    binomial(n_sites_ - targetnup, targetndn);

    // Try allocating upspin / downspin vectors
    std::vector<bit_t> baseupspins;
    std::vector<bit_t> basednspins;
    std::vector<bit_t> targetupspins;
    std::vector<bit_t> targetdnspins;
    try {
      baseupspins.resize(dim_);
      basednspins.resize(dim_);
      targetupspins.resize(targetdim);
      targetdnspins.resize(targetdim);
    } catch (...) {
      std::cerr << "Error: Could not allocate upspin/downspin "
                << "vectors for TJModel!" << std::endl
                << std::flush;
      exit(EXIT_FAILURE);
    }

    // Fill vectors holding upspin/downspin configurations
    int64 idx = 0;
    auto hs_baseupspins = BasisSpinHalf<bit_t>(n_sites_, basenup);
    auto hs_holes_in_baseups =
        BasisSpinHalf<bit_t>(n_sites_ - basenup, basendn);
    for (auto ups : hs_baseupspins)
      for (auto holes : hs_holes_in_baseups) {
        auto dns = up_hole_to_down(ups.spins, holes.spins);
        baseupspins[idx] = ups.spins;
        basednspins[idx] = dns;
        ++idx;
      }
    assert(idx == dim_);

    idx = 0;
    auto hs_targetupspins = BasisSpinHalf<bit_t>(n_sites_, targetnup);
    auto hs_holes_in_targetups =
        BasisSpinHalf<bit_t>(n_sites_ - targetnup, targetndn);
    for (auto ups : hs_targetupspins)
      for (auto holes : hs_holes_in_targetups) {
        auto dns = up_hole_to_down(ups.spins, holes.spins);
        targetupspins[idx] = ups.spins;
        targetdnspins[idx] = dns;
        ++idx;
      }
    assert(idx == targetdim);

    // Try allocating the matrix
    lila::Matrix<double> sMinus;
    try {
      sMinus.resize(targetdim, dim_);
      Zeros(sMinus);
    } catch (...) {
      std::cerr << "Error: Could not allocate matrix for TJModel!" << std::endl
                << std::flush;
      exit(EXIT_FAILURE);
    }

    // Define lambda function to get indices for target space
    auto index_of_up_dn = [&targetupspins, &targetdnspins](bit_t const &ups,
                                                           bit_t const &dns) {
      // Binary search the new indices
      auto up_bounds =
          std::equal_range(targetupspins.begin(), targetupspins.end(), ups);
      auto up_begin = std::distance(targetupspins.begin(), up_bounds.first);
      auto up_end = std::distance(targetupspins.begin(), up_bounds.second);
      auto it = std::lower_bound(targetdnspins.begin() + up_begin,
                                 targetdnspins.begin() + up_end, dns);
      return std::distance(targetdnspins.begin(), it);
    };

    // Iterate through spin configurations
    bit_t flipmask = ((bit_t)1 << siteIndex);

    for (int64 idx = 0; idx < dim_; ++idx) {
      bit_t ups = baseupspins[idx];
      bit_t dns = basednspins[idx];

      bool up = gbit(ups, siteIndex);

      if (up) {
        bit_t flipped_ups = ups ^ flipmask;
        bit_t flipped_dns = dns ^ flipmask;
        int64 target_idx = index_of_up_dn(flipped_ups, flipped_dns);
        sMinus(target_idx, idx) = 1;
      }
    }

    return sMinus;

  } else {

    lila::Matrix<double> sMinus;
    sMinus.resize(1, 1);
    Zeros(sMinus);
    return sMinus;
  }
}

template class TJModel<double, uint32, uint64>;
template class TJModel<double, uint64, uint64>;

template class TJModel<complex, uint32, uint64>;
template class TJModel<complex, uint64, uint64>;

} // namespace hydra
