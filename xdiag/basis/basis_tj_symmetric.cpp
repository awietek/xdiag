// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "basis_tj_symmetric.hpp"

#include <cmath>

#include <xdiag/bits/bitset.hpp>
#include <xdiag/bits/zero_one.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag::basis {

// Coupled fermionic norm of the symmetric tJ state (ups, dns): identical to the
// electron case (the no-double-occupancy constraint only restricts WHICH (ups,
// dns) are stored, not the norm formula). The up-stabilizer `stab` fixes `ups`,
// so only its elements that also fix `dns` contribute, each weighted by the
// irrep character and the fermi sign fermi_up XOR fermi_dn.
template <typename enumeration_t>
static double
coupled_norm(typename enumeration_t::bit_t ups,
             typename enumeration_t::bit_t dns,
             symmetries::SitePermutation const &action,
             std::vector<int64_t> const &stab, arma::cx_vec const &characters,
             symmetries::FermiTable<enumeration_t> const &fermi_up,
             symmetries::FermiTable<enumeration_t> const &fermi_dn) {
  complex amplitude = 0.0;
  for (int64_t sym : stab) {
    if (action.apply(sym, dns) == dns) {
      bool fu = fermi_up.sign(sym, ups);
      bool fd = fermi_dn.sign(sym, dns);
      amplitude += (fu == fd) ? characters(sym) : -characters(sym);
    }
  }
  return std::sqrt(std::abs(amplitude));
}

template <typename enumeration_t>
BasistJSymmetric<enumeration_t>::BasistJSymmetric(
    enumeration_t const &enum_up, enumeration_t const &enum_dn,
    PermutationGroup const &group, Vector const &characters) try
    : basis_up_(enum_up, group,
                Vector(arma::vec(group.size(), arma::fill::ones)),
                /*fermionic=*/false),
      basis_dn_(enum_dn), characters_(characters), action_(group),
      fermi_up_(enum_up, group), fermi_dn_(enum_dn, group) {

  if (enum_up.n() != enum_dn.n()) {
    XDIAG_THROW(fmt::format(
        "up and dn basis have different number of sites (up: {}), (dn: {})",
        enum_up.n(), enum_dn.n()));
  }
  if (enum_up.n() != group.nsites()) {
    XDIAG_THROW("nsites of the enumeration does not match the nsites of the "
                "PermutationGroup");
  }

  int64_t nsites = enum_up.n();
  bit_t sitesmask{};

  // Compressed-dn rank over Combinations(nsites - nup, ndn). Only built for the
  // number-conserving integral case; the branch is discarded by if constexpr
  // otherwise so enum.k() / LinTable<BitsetDynamic> are never instantiated for
  // Subsets / BitsetDynamic.
  if constexpr (use_compressed_index_) {
    sitesmask = (bit_t)((((bit_t)1) << nsites) - 1);
    if constexpr (!combinatorics::is_subsets_v<enumeration_t>) {
      int64_t nup = enum_up.k();
      int64_t ndn = enum_dn.k();
      lintable_dnsc_ = combinatorics::LinTable<bit_t>(nsites - nup, ndn);
    }
  }

  int64_t n_rep_ups = basis_up_.size();
  arma::cx_vec chars = characters_.as<arma::cx_vec>();

  if constexpr (use_compressed_index_) {
    extractors_.reserve(n_rep_ups);
  }
  ups_offset_.resize(n_rep_ups);
  // (start, length) of each up representative's dn block in the flat storage;
  // converted to spans below, once the storage has stopped reallocating.
  std::vector<std::pair<span_size_t, span_size_t>> dns_limits(n_rep_ups);

  // Unlike BasisElectronSymmetric there is NO shared front block: the
  // no-double-occupancy constraint `ups & dns == 0` couples the allowed dn
  // configurations to the up representative, so EVERY up representative
  // materialises its own dn block (filtered by the constraint). One-time
  // construction, kept serial.
  size_ = 0;
  for (int64_t idx_up = 0; idx_up < n_rep_ups; ++idx_up) {
    bit_t ups = basis_up_[idx_up];
    ups_offset_[idx_up] = size_;
    span_size_t start = dns_storage_.size();

    // precompile this up rep's compression mask (its non-up sites)
    if constexpr (use_compressed_index_) {
      extractors_.emplace_back((bit_t)((~ups) & sitesmask));
    }

    if (stab_size(idx_up) == 1) {
      // trivial up-stabilizer: every dn with no double occupancy, norm 1. The
      // dn enumeration is ascending, so the block stays ascending.
      for (bit_t dns : basis_dn_) {
        if (bits::iszero(dns & ups)) {
          dns_storage_.push_back(dns);
          norms_storage_.push_back(1.0);
        }
      }
    } else {
      // non-trivial up-stabilizer: materialise the dn-rep block, dropping
      // double-occupied configurations and zero-norm representatives. The
      // stabilizer of the representative equals syms_ups(ups) since ups is a
      // rep here.
      std::vector<int64_t> stab = syms_ups(ups);
      for (bit_t dns : basis_dn_) {
        if (!bits::iszero(dns & ups)) {
          continue; // no double occupancy
        }
        bit_t dns_rep = symmetries::representative_subset(dns, action_, stab);
        if (dns == dns_rep) {
          double nrm = coupled_norm<enumeration_t>(ups, dns, action_, stab,
                                                   chars, fermi_up_, fermi_dn_);
          if (nrm > 1e-6) {
            dns_storage_.push_back(dns);
            norms_storage_.push_back(nrm);
          }
        }
      }
    }
    span_size_t length = dns_storage_.size() - start;
    dns_limits[idx_up] = {start, length};
    size_ += (int64_t)length;
  }

  // Storage is now final: materialise the per-rep spans (cannot dangle now).
  dns_for_ups_rep_.reserve(n_rep_ups);
  norms_for_ups_rep_.reserve(n_rep_ups);
  for (int64_t idx_up = 0; idx_up < n_rep_ups; ++idx_up) {
    auto [start, length] = dns_limits[idx_up];
    dns_for_ups_rep_.push_back({dns_storage_.data() + start, length});
    norms_for_ups_rep_.push_back({norms_storage_.data() + start, length});
  }
}
XDIAG_CATCH

template <typename enumeration_t>
int64_t BasistJSymmetric<enumeration_t>::size() const {
  return size_;
}

template <typename enumeration_t>
int64_t BasistJSymmetric<enumeration_t>::nsites() const {
  return basis_dn_.nsites();
}

template <typename enumeration_t>
int64_t BasistJSymmetric<enumeration_t>::index(bit_t ups, bit_t dns) const {
  auto [idx, sym, norm_out, fermi] = representative_data_fermi(ups, dns);
  (void)sym;
  (void)norm_out;
  (void)fermi;
  return idx;
}

template <typename enumeration_t>
int64_t
BasistJSymmetric<enumeration_t>::index(ProductState const &pstate) const {
  auto [ups, dns] = pair_from_pstate<bit_t>(pstate, nsites());
  return index(ups, dns);
}

template <typename enumeration_t>
std::vector<int64_t>
BasistJSymmetric<enumeration_t>::syms_ups(bit_t ups) const {
  auto [raw, s0, nrm] = basis_up_.representative_data(ups);
  (void)s0;
  (void)nrm;
  bit_t rep = basis_up_[raw - 1];
  std::vector<int64_t> syms;
  for (int64_t g = 0; g < action_.size(); ++g) {
    if (action_.apply(g, ups) == rep) {
      syms.push_back(g);
    }
  }
  return syms;
}

template <typename enumeration_t>
std::tuple<int64_t, bool, int64_t>
BasistJSymmetric<enumeration_t>::index_dns_fermi_sym(
    bit_t dns, std::vector<int64_t> const &syms,
    gsl::span<bit_t const> dnss) const {
  auto [rep_dns, rep_sym] =
      symmetries::representative_sym_subset(dns, action_, syms);
  auto it = std::lower_bound(dnss.begin(), dnss.end(), rep_dns);
  if ((it != dnss.end()) && (*it == rep_dns)) {
    return {(int64_t)std::distance(dnss.begin(), it),
            fermi_bool_dns(rep_sym, dns), rep_sym};
  } else {
    return {-1, false, rep_sym};
  }
}

template <typename enumeration_t>
std::tuple<int64_t, int64_t, double, bool>
BasistJSymmetric<enumeration_t>::representative_data_fermi(bit_t ups,
                                                          bit_t dns) const {
  auto [raw, s0, nrm] = basis_up_.representative_data(ups);
  (void)nrm;
  int64_t idx_ups = raw - 1; // with all-ones characters raw is always >= 1
  if (stab_size(idx_ups) == 1) {
    auto [idx_dn, fermi_dn] =
        index_dns_fermi(dns, s0, idx_ups, dns_for_ups_rep(idx_ups));
    if (idx_dn < 0) {
      return {-1, s0, 0.0, false};
    }
    bool fermi = fermi_bool_ups(s0, ups) ^ fermi_dn;
    return {ups_offset_[idx_ups] + idx_dn, s0, 1.0, fermi};
  } else {
    std::vector<int64_t> syms = syms_ups(ups);
    auto [idx_dn, fermi_dn, sym] =
        index_dns_fermi_sym(dns, syms, dns_for_ups_rep(idx_ups));
    if (idx_dn < 0) {
      return {-1, sym, 0.0, false};
    }
    double norm_out = norms_for_ups_rep(idx_ups)[idx_dn];
    bool fermi = fermi_bool_ups(sym, ups) ^ fermi_dn;
    return {ups_offset_[idx_ups] + idx_dn, sym, norm_out, fermi};
  }
}

template <typename enumeration_t>
typename BasistJSymmetric<enumeration_t>::iterator_t
BasistJSymmetric<enumeration_t>::begin() const {
  return iterator_t(*this, true);
}

template <typename enumeration_t>
typename BasistJSymmetric<enumeration_t>::iterator_t
BasistJSymmetric<enumeration_t>::end() const {
  return iterator_t(*this, false);
}

template <typename enumeration_t>
bool BasistJSymmetric<enumeration_t>::operator==(
    BasistJSymmetric<enumeration_t> const &rhs) const {
  return (basis_up_ == rhs.basis_up_) && (basis_dn_ == rhs.basis_dn_) &&
         (group() == rhs.group()) && (characters_ == rhs.characters_);
}

template <typename enumeration_t>
bool BasistJSymmetric<enumeration_t>::operator!=(
    BasistJSymmetric<enumeration_t> const &rhs) const {
  return !operator==(rhs);
}

template <typename enumeration_t>
BasistJSymmetricIterator<enumeration_t>::BasistJSymmetricIterator(
    BasistJSymmetric<enumeration_t> const &basis, bool begin)
    : basis_(&basis), up_idx_(begin ? 0 : basis.basis_up().size()),
      dn_idx_(0) {
  if ((basis.basis_up().size() > 0) && begin) {
    dnss_ = basis.dns_for_ups_rep(0);
    // skip leading up representatives whose dn block is empty
    while ((up_idx_ < basis.basis_up().size()) && (dnss_.size() == 0)) {
      ++up_idx_;
      if (up_idx_ < basis.basis_up().size()) {
        dnss_ = basis.dns_for_ups_rep(up_idx_);
      }
    }
  }
}

template <typename enumeration_t>
BasistJSymmetricIterator<enumeration_t> &
BasistJSymmetricIterator<enumeration_t>::operator++() {
  ++dn_idx_;
  if (dn_idx_ == (int64_t)dnss_.size()) {
    dn_idx_ = 0;
    do {
      ++up_idx_;
      if (up_idx_ == basis_->basis_up().size()) {
        return *this;
      }
      dnss_ = basis_->dns_for_ups_rep(up_idx_);
    } while (dnss_.size() == 0);
  }
  return *this;
}

template <typename enumeration_t>
std::pair<typename BasistJSymmetricIterator<enumeration_t>::bit_t,
          typename BasistJSymmetricIterator<enumeration_t>::bit_t>
BasistJSymmetricIterator<enumeration_t>::operator*() const {
  return {basis_->basis_up()[up_idx_], dnss_[dn_idx_]};
}

template <typename enumeration_t>
bool BasistJSymmetricIterator<enumeration_t>::operator!=(
    BasistJSymmetricIterator<enumeration_t> const &rhs) const {
  return (up_idx_ != rhs.up_idx_) || (dn_idx_ != rhs.dn_idx_);
}

} // namespace xdiag::basis

using namespace xdiag::combinatorics;
using namespace xdiag::bits;

#define INSTANTIATE_BASIS_TJ_SYMMETRIC(ENUMERATION)                            \
  template class xdiag::basis::BasistJSymmetric<ENUMERATION>;                  \
  template class xdiag::basis::BasistJSymmetricIterator<ENUMERATION>;

INSTANTIATE_BASIS_TJ_SYMMETRIC(Subsets<uint32_t>)
INSTANTIATE_BASIS_TJ_SYMMETRIC(Subsets<uint64_t>)
INSTANTIATE_BASIS_TJ_SYMMETRIC(LinTable<uint32_t>)
INSTANTIATE_BASIS_TJ_SYMMETRIC(LinTable<uint64_t>)
INSTANTIATE_BASIS_TJ_SYMMETRIC(Combinations<uint32_t>)
INSTANTIATE_BASIS_TJ_SYMMETRIC(Combinations<uint64_t>)
INSTANTIATE_BASIS_TJ_SYMMETRIC(Combinations<BitsetDynamic>)

#undef INSTANTIATE_BASIS_TJ_SYMMETRIC
