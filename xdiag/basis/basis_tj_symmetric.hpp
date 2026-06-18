// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include <xdiag/extern/gsl/span>

#include <type_traits>

#include <xdiag/basis/basis.hpp>
#include <xdiag/basis/basis_onthefly.hpp>
#include <xdiag/basis/basis_symmetric.hpp>
#include <xdiag/bits/mask_compressor.hpp>
#include <xdiag/bits/zero_one.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/math/vector.hpp>
#include <xdiag/symmetries/action/representative.hpp>
#include <xdiag/symmetries/action/site_permutation.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/tables/fermi_table.hpp>
#include <xdiag/utils/type_name.hpp>

namespace xdiag::basis {

template <typename enumeration_tt> class BasistJSymmetricIterator;

// tJ basis (local dimension d = 3: empty / up / dn, no double occupancy) with
// permutation symmetry. Same coupled scheme as BasisElectronSymmetric -- the up
// sector is reduced to orbit representatives under the FULL group (bosonically,
// all-ones characters, so the stored up norm is sqrt(|up-stabilizer|)), and for
// each up representative the dn configurations are symmetrised under that
// representative's up-stabilizer with the coupled fermionic norm. The symmetric
// state index is ups_offset(idx_ups) + idx_dns_within_block.
//
// The ONE structural difference from BasisElectronSymmetric is the dn storage.
// Electron shares a single front block of ALL dn states for every up
// representative with a trivial up-stabilizer (double occupancy is allowed, so
// the dn block is up-independent). The tJ no-double-occupancy constraint
// `ups & dns == 0` couples the allowed dn configurations to the up
// representative, so there is NO shared front: EVERY up representative (trivial
// and non-trivial up-stabilizer) materialises its own dn block, filtered by
// `(ups_rep & dns) == 0`. (That filter is complete: a site permutation g acts on
// both sectors, so (g.ups) & (g.dns) = g.(ups & dns); the constraint is
// permutation-invariant, hence checking it on the representative suffices.)
//
// The dn lookups (index_dns, index_dns_fermi) return -1 when the dn is absent
// from the up-rep block -- which is exactly how a double-occupancy intermediate
// is dropped. The symmetric matrix kernels are shared with the electron block
// (templated on the basis type); they only differ from electron in that these
// lookups can return -1 (electron's never do), and the kernels already guard on
// idx >= 0.
//
// Fast dn index: although the block is up-rep-specific, the compressed dn (dns
// packed into the nsites - nup non-up sites) is order-preserving, and so is the
// dn enumeration order the block is built in, so the block position is a direct
// function of the compressed dn. The compression mask `not_ups` is constant
// across a whole inner dn loop, so it is precompiled once per up rep into a
// bits::MaskCompressor (branchless O(log wordsize), no BMI2 needed). Then:
//   - no number conservation (Subsets): the compressed value IS the block index;
//   - number-conserving (Combinations / LinTable): one O(1) LinTable rank of the
//     compressed dn over Combinations(nsites - nup, ndn).
// This covers the integral-bit cases (uint32/uint64); the large BitsetDynamic
// case (no feasible table / integer compress) falls back to a binary search of
// the stored block.
template <typename enumeration_tt>
class BasistJSymmetric : public BasisType<BasistJSymmetric<enumeration_tt>> {
public:
  using enumeration_t = enumeration_tt;
  using bit_t = typename enumeration_t::bit_t;
  using iterator_t = BasistJSymmetricIterator<enumeration_t>;
  using span_size_t = gsl::span<int64_t const>::size_type;
  static constexpr std::string_view type_name =
      utils::get_type_name<BasistJSymmetric<enumeration_t>>();

  BasistJSymmetric() = default;
  BasistJSymmetric(enumeration_t const &enum_up, enumeration_t const &enum_dn,
                   PermutationGroup const &group, Vector const &characters);

  int64_t size() const override;
  int64_t nsites() const override;
  constexpr int64_t d() const { return 3; }

  // Linear index of the symmetric basis state reached by symmetrising the raw
  // configuration (ups, dns); -1 if it has zero norm or is not in the basis
  // (which includes a double-occupied configuration).
  int64_t index(bit_t ups, bit_t dns) const;

  iterator_t begin() const;
  iterator_t end() const;

  // The up representatives (BasisSymmetric) and dn enumeration (BasisOnTheFly).
  BasisSymmetric<enumeration_t> const &basis_up() const { return basis_up_; }
  BasisOnTheFly<enumeration_t> const &basis_dn() const { return basis_dn_; }

  PermutationGroup const &group() const { return action_.group(); }
  Vector const &characters() const { return characters_; }

  bool operator==(BasistJSymmetric<enumeration_t> const &rhs) const;
  bool operator!=(BasistJSymmetric<enumeration_t> const &rhs) const;

  // ---- coupled data used by the symmetric matrix kernels ----------------

  inline int64_t ups_offset(int64_t idx_ups) const {
    return ups_offset_[idx_ups];
  }

  // Order of the up-stabilizer of representative idx_ups (== 1 for a trivial
  // stabilizer). With all-ones characters in basis_up_ the stored up norm is
  // sqrt(|stabilizer|), so this is O(1).
  inline int64_t stab_size(int64_t idx_ups) const {
    double n = basis_up_.norm(idx_ups);
    return (int64_t)std::llround(n * n);
  }

  // All symmetries mapping `ups` to its representative (the up-stabilizer coset).
  std::vector<int64_t> syms_ups(bit_t ups) const;

  inline gsl::span<bit_t const> dns_for_ups_rep(int64_t idx_ups) const {
    return dns_for_ups_rep_[idx_ups];
  }
  inline gsl::span<double const> norms_for_ups_rep(int64_t idx_ups) const {
    return norms_for_ups_rep_[idx_ups];
  }

  // O(1) fermi signs of a single symmetry acting on an up / dn configuration.
  inline bool fermi_bool_ups(int64_t sym, bit_t ups) const {
    return fermi_up_.sign(sym, ups);
  }
  inline bool fermi_bool_dns(int64_t sym, bit_t dns) const {
    return fermi_dn_.sign(sym, dns);
  }

  // Whether the precompiled-compress dn index applies: the integral-bit cases
  // (Combinations / LinTable / Subsets over uint*). BitsetDynamic (no feasible
  // table / integer compress) falls back to a binary search of the stored block.
  static constexpr bool use_compressed_index_ = std::is_integral_v<bit_t>;

  // Block index of `dns` in the dn block of up representative `idx_up`; -1 if
  // absent (which is how a double-occupied dn output is dropped). For the
  // integral case this is the precompiled-compress lookup (no number
  // conservation: the compressed value is the index; number-conserving: its
  // LinTable rank); otherwise a binary search of the ascending block `dnss`.
  inline int64_t index_dns(bit_t dns, int64_t idx_up,
                           gsl::span<bit_t const> dnss) const {
    if constexpr (use_compressed_index_) {
      (void)dnss;
      bits::MaskCompressor<bit_t> const &ex = extractors_[idx_up];
      if ((dns & ex.mask()) != dns) {
        return -1; // a bit outside the non-up sites -> double occupancy
      }
      bit_t c = ex(dns);
      if constexpr (combinatorics::is_subsets_v<enumeration_t>) {
        return (int64_t)c; // compressed value is the block index
      } else {
        return lintable_dnsc_.index(c); // combinatorial rank in the block
      }
    } else {
      (void)idx_up;
      auto it = std::lower_bound(dnss.begin(), dnss.end(), dns);
      if ((it != dnss.end()) && (*it == dns)) {
        return (int64_t)std::distance(dnss.begin(), it);
      }
      return -1;
    }
  }

  // dn index (within up representative `idx_up`'s block) and fermi sign of `dns`
  // mapped by the single symmetry `sym` (the trivial up-stabilizer path).
  // idx == -1 if absent.
  inline std::pair<int64_t, bool>
  index_dns_fermi(bit_t dns, int64_t sym, int64_t idx_up,
                  gsl::span<bit_t const> dnss) const {
    bit_t dns_rep = action_.apply(sym, dns);
    return {index_dns(dns_rep, idx_up, dnss), fermi_bool_dns(sym, dns)};
  }

  // dn index, fermi sign and chosen symmetry for the non-trivial up-stabilizer.
  // `syms` is the up-stabilizer coset; `dnss` is the (ascending) dn-rep block.
  // idx_dn == -1 means `dns` maps to an absent dn representative.
  std::tuple<int64_t, bool, int64_t>
  index_dns_fermi_sym(bit_t dns, std::vector<int64_t> const &syms,
                      gsl::span<bit_t const> dnss) const;

  // Unified "symmetrise an arbitrary (ups, dns)" primitive. Returns
  // {index, sym, norm_out, fermi_total}; index == -1 means the state has zero
  // norm or is not in the basis (includes double occupancy).
  std::tuple<int64_t, int64_t, double, bool>
  representative_data_fermi(bit_t ups, bit_t dns) const;

private:
  // up sector reduced to orbit representatives (trivial all-ones characters)
  BasisSymmetric<enumeration_t> basis_up_;
  // full dn enumeration (iterated once per up rep to build that rep's dn block)
  BasisOnTheFly<enumeration_t> basis_dn_;

  Vector characters_;                  // the actual 1-D irrep characters
  symmetries::SitePermutation action_; // owns a copy of the group

  symmetries::FermiTable<enumeration_t> fermi_up_;
  symmetries::FermiTable<enumeration_t> fermi_dn_;

  // Per-up-rep precompiled compressors (mask = that rep's non-up sites). Built
  // and used only when use_compressed_index_; empty otherwise.
  std::vector<bits::MaskCompressor<bit_t>> extractors_;
  // Compressed-dn rank over Combinations(nsites - nup, ndn). Built and used only
  // for the number-conserving integral case; default-constructed otherwise.
  combinatorics::LinTable<bit_t> lintable_dnsc_;

  // dn-rep blocks: one bespoke block per up representative (no shared front, as
  // the allowed dns depend on the up rep). The per-rep spans view into the flat
  // storage, materialised once after the storage is finalised.
  std::vector<bit_t> dns_storage_;
  std::vector<double> norms_storage_;
  std::vector<gsl::span<bit_t const>> dns_for_ups_rep_;
  std::vector<gsl::span<double const>> norms_for_ups_rep_;
  std::vector<int64_t> ups_offset_;

  int64_t size_ = 0;
};

// Iterates (ups_rep, dns_rep) pairs: up representatives outer, dn block inner.
template <typename enumeration_tt> class BasistJSymmetricIterator {
public:
  using enumeration_t = enumeration_tt;
  using bit_t = typename enumeration_t::bit_t;

  BasistJSymmetricIterator() = default;
  BasistJSymmetricIterator(BasistJSymmetric<enumeration_t> const &basis,
                           bool begin);
  BasistJSymmetricIterator &operator++();
  std::pair<bit_t, bit_t> operator*() const;
  bool operator!=(BasistJSymmetricIterator<enumeration_t> const &rhs) const;

private:
  BasistJSymmetric<enumeration_t> const *basis_ = nullptr;
  int64_t up_idx_ = 0;
  int64_t dn_idx_ = 0;
  gsl::span<bit_t const> dnss_;
};

} // namespace xdiag::basis
