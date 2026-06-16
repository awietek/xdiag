// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "basis_electron_symmetric.hpp"

#include <cmath>

#include <xdiag/bits/bitset.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag::basis {

// Coupled fermionic norm of the symmetric electron state (ups, dns): the
// stabilizer subset `stab` already fixes `ups`, so only its elements that also
// fix `dns` (the joint stabilizer) contribute, each weighted by the irrep
// character and the fermi sign fermi_up XOR fermi_dn, read from the precomputed
// fermi tables. Complex characters are handled uniformly (a real irrep embeds
// as a real-valued complex vector).
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
BasisElectronSymmetric<enumeration_t>::BasisElectronSymmetric(
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

  int64_t n_rep_ups = basis_up_.size();
  int64_t size_dn = basis_dn_.size();
  arma::cx_vec chars = characters_.as<arma::cx_vec>();

  ups_offset_.resize(n_rep_ups);
  // (start, length) of each up representative's dn block in the flat storage;
  // converted to spans below, once the storage has stopped reallocating.
  std::vector<std::pair<span_size_t, span_size_t>> dns_limits(n_rep_ups);

  // Shared front block: all dn states (norm 1), used by every up representative
  // with a trivial up-stabilizer. The dn enumeration is ascending, so the front
  // doubles as the (ascending) search block for the trivial case.
  for (bit_t dns : basis_dn_) {
    dns_storage_.push_back(dns);
    norms_storage_.push_back(1.0);
  }

  // One-time construction, kept serial: a trivial up-stabilizer (the common
  // case) is O(1) -- it just points at the shared front -- and only the rare
  // non-trivial up-stabilizers do the O(size_dn) dn scan. It parallelises over
  // representatives if it ever becomes a bottleneck (each rep fills a disjoint
  // block), but that is not worth the added machinery here.
  size_ = 0;
  for (int64_t idx_up = 0; idx_up < n_rep_ups; ++idx_up) {
    bit_t ups = basis_up_[idx_up];
    ups_offset_[idx_up] = size_;

    if (stab_size(idx_up) == 1) {
      // trivial up-stabilizer: point at the shared front block
      dns_limits[idx_up] = {(span_size_t)0, (span_size_t)size_dn};
      size_ += size_dn;
    } else {
      // non-trivial up-stabilizer: materialise the dn-rep block. The stabilizer
      // of the representative equals syms_ups(ups) since ups is a rep here.
      std::vector<int64_t> stab = syms_ups(ups);
      span_size_t start = dns_storage_.size();
      for (bit_t dns : basis_dn_) {
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
      span_size_t length = dns_storage_.size() - start;
      dns_limits[idx_up] = {start, length};
      size_ += (int64_t)length;
    }
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
int64_t BasisElectronSymmetric<enumeration_t>::size() const {
  return size_;
}

template <typename enumeration_t>
int64_t BasisElectronSymmetric<enumeration_t>::nsites() const {
  return basis_dn_.nsites();
}

template <typename enumeration_t>
int64_t BasisElectronSymmetric<enumeration_t>::index(bit_t ups,
                                                     bit_t dns) const {
  auto [idx, sym, norm_out, fermi] = representative_data_fermi(ups, dns);
  (void)sym;
  (void)norm_out;
  (void)fermi;
  return idx;
}

template <typename enumeration_t>
std::vector<int64_t>
BasisElectronSymmetric<enumeration_t>::syms_ups(bit_t ups) const {
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
BasisElectronSymmetric<enumeration_t>::index_dns_fermi_sym(
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
BasisElectronSymmetric<enumeration_t>::representative_data_fermi(
    bit_t ups, bit_t dns) const {
  auto [raw, s0, nrm] = basis_up_.representative_data(ups);
  (void)nrm;
  int64_t idx_ups = raw - 1; // with all-ones characters raw is always >= 1
  if (stab_size(idx_ups) == 1) {
    auto [idx_dn, fermi_dn] = index_dns_fermi(dns, s0);
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
typename BasisElectronSymmetric<enumeration_t>::iterator_t
BasisElectronSymmetric<enumeration_t>::begin() const {
  return iterator_t(*this, true);
}

template <typename enumeration_t>
typename BasisElectronSymmetric<enumeration_t>::iterator_t
BasisElectronSymmetric<enumeration_t>::end() const {
  return iterator_t(*this, false);
}

template <typename enumeration_t>
bool BasisElectronSymmetric<enumeration_t>::operator==(
    BasisElectronSymmetric<enumeration_t> const &rhs) const {
  return (basis_up_ == rhs.basis_up_) && (basis_dn_ == rhs.basis_dn_) &&
         (group() == rhs.group()) && (characters_ == rhs.characters_);
}

template <typename enumeration_t>
bool BasisElectronSymmetric<enumeration_t>::operator!=(
    BasisElectronSymmetric<enumeration_t> const &rhs) const {
  return !operator==(rhs);
}

template <typename enumeration_t>
BasisElectronSymmetricIterator<enumeration_t>::BasisElectronSymmetricIterator(
    BasisElectronSymmetric<enumeration_t> const &basis, bool begin)
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
BasisElectronSymmetricIterator<enumeration_t> &
BasisElectronSymmetricIterator<enumeration_t>::operator++() {
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
std::pair<typename BasisElectronSymmetricIterator<enumeration_t>::bit_t,
          typename BasisElectronSymmetricIterator<enumeration_t>::bit_t>
BasisElectronSymmetricIterator<enumeration_t>::operator*() const {
  return {basis_->basis_up()[up_idx_], dnss_[dn_idx_]};
}

template <typename enumeration_t>
bool BasisElectronSymmetricIterator<enumeration_t>::operator!=(
    BasisElectronSymmetricIterator<enumeration_t> const &rhs) const {
  return (up_idx_ != rhs.up_idx_) || (dn_idx_ != rhs.dn_idx_);
}

} // namespace xdiag::basis

using namespace xdiag::combinatorics;
using namespace xdiag::bits;

#define INSTANTIATE_BASIS_ELECTRON_SYMMETRIC(ENUMERATION)                      \
  template class xdiag::basis::BasisElectronSymmetric<ENUMERATION>;            \
  template class xdiag::basis::BasisElectronSymmetricIterator<ENUMERATION>;

INSTANTIATE_BASIS_ELECTRON_SYMMETRIC(Subsets<uint32_t>)
INSTANTIATE_BASIS_ELECTRON_SYMMETRIC(Subsets<uint64_t>)
INSTANTIATE_BASIS_ELECTRON_SYMMETRIC(LinTable<uint32_t>)
INSTANTIATE_BASIS_ELECTRON_SYMMETRIC(LinTable<uint64_t>)
INSTANTIATE_BASIS_ELECTRON_SYMMETRIC(Combinations<uint32_t>)
INSTANTIATE_BASIS_ELECTRON_SYMMETRIC(Combinations<uint64_t>)
INSTANTIATE_BASIS_ELECTRON_SYMMETRIC(Combinations<BitsetDynamic>)

#undef INSTANTIATE_BASIS_ELECTRON_SYMMETRIC
