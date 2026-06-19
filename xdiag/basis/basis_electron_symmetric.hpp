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

#include <extern/gsl/span>

#include <xdiag/basis/basis.hpp>
#include <xdiag/basis/basis_onthefly.hpp>
#include <xdiag/basis/basis_symmetric.hpp>
#include <xdiag/math/vector.hpp>
#include <xdiag/symmetries/action/representative.hpp>
#include <xdiag/symmetries/action/site_permutation.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/tables/fermi_table.hpp>
#include <xdiag/utils/type_name.hpp>

namespace xdiag::basis {

template <typename enumeration_tt> class BasisElectronSymmetricIterator;

// Spinful electron basis (local dimension d = 4) with permutation symmetry,
// templated on the enumeration type exactly like BasisElectron: Subsets ⇒ no
// number conservation, Combinations / LinTable ⇒ fixed (nup, ndn).
//
// The group acts on (ups, dns) jointly, so the two sectors cannot be
// symmetrised independently. A single representative table over the joint
// (ups, dns) space would cost O(D_full) memory, so we use the memory-efficient
// coupled scheme:
//   - the up sector is reduced to its orbit representatives under the FULL
//     group. This reuses BasisSymmetric (basis_up()) with TRIVIAL all-ones
//     characters, which drops nothing: the *fermionic* up norm of a stabilizer
//     is |sum_{g in Stab} eps_up(g)|, an orthogonality sum that vanishes
//     whenever the up-stabilizer's fermi character eps_up is non-trivial -- yet
//     such up orbits still host valid electron states (a dn that breaks the
//     stabilizer gives a non-zero coupled norm). So the ups must be reduced
//     bosonically, keeping every orbit's representative as a label. A bonus of
//     the all-ones characters: the stored up norm equals sqrt(|up-stabilizer|),
//     which is how stab_size() recovers the stabilizer order in O(1).
//   - for each up representative the dn configurations are symmetrised under
//     that representative's up-stabilizer subgroup, with the coupled fermionic
//     norm (joint stabilizer, weighted by fermi_up XOR fermi_dn).
//   - a symmetric state index is ups_offset(idx_ups) + idx_dns_within_block.
//
// As a memory optimisation up representatives with a TRIVIAL up-stabilizer
// (the common case) share a single front block of all dn states (norm 1, the
// full dn enumeration basis_dn(), stored once); only non-trivial up-stabilizers
// materialise a bespoke dn-rep block.
//
// Because the ups table is bosonic it does not store the ups fermi sign (the
// representative table only does so for fermionic tables). The ups fermi sign
// is therefore recomputed from the chosen symmetry; in the kernels this is
// hoisted out of the inner dn loop (once per up representative), so its cost is
// negligible next to the per-element dn work.
//
// Operations on the up representatives (count, representative, representative
// index/symmetry) and the dn enumeration (index) are NOT re-wrapped here: the
// kernels reach through basis_up() (a BasisSymmetric) and basis_dn() (a
// BasisOnTheFly), just as the non-symmetric kernels reach through
// BasisElectron::basis_up()/basis_dn(). Only the genuinely coupled data
// (offsets, dn blocks, up-stabilizers, joint symmetrisation) lives here.
template <typename enumeration_tt>
class BasisElectronSymmetric
    : public BasisType<BasisElectronSymmetric<enumeration_tt>> {
public:
  using enumeration_t = enumeration_tt;
  using bit_t = typename enumeration_t::bit_t;
  using iterator_t = BasisElectronSymmetricIterator<enumeration_t>;
  using span_size_t = gsl::span<int64_t const>::size_type;
  static constexpr std::string_view type_name =
      utils::get_type_name<BasisElectronSymmetric<enumeration_t>>();

  BasisElectronSymmetric() = default;
  BasisElectronSymmetric(enumeration_t const &enum_up,
                         enumeration_t const &enum_dn,
                         PermutationGroup const &group,
                         Vector const &characters);

  int64_t size() const override;
  int64_t nsites() const override;
  constexpr int64_t d() const { return 4; }

  // Linear index of the symmetric basis state reached by symmetrising the raw
  // configuration (ups, dns); -1 if it has zero norm (not in the basis).
  int64_t index(bit_t ups, bit_t dns) const;

  iterator_t begin() const;
  iterator_t end() const;

  // The up representatives (BasisSymmetric) and dn enumeration (BasisOnTheFly).
  // Source of n_rep_ups (size()), rep_ups (operator[]), index_ups / rep_sym
  // (representative_data) and index_dns (index) -- reached through directly.
  BasisSymmetric<enumeration_t> const &basis_up() const { return basis_up_; }
  BasisOnTheFly<enumeration_t> const &basis_dn() const { return basis_dn_; }

  PermutationGroup const &group() const { return action_.group(); }
  Vector const &characters() const { return characters_; }

  bool operator==(BasisElectronSymmetric<enumeration_t> const &rhs) const;
  bool operator!=(BasisElectronSymmetric<enumeration_t> const &rhs) const;

  // ---- coupled data used by the symmetric matrix kernels ----------------

  inline int64_t ups_offset(int64_t idx_ups) const {
    return ups_offset_[idx_ups];
  }

  // Order of the up-stabilizer of representative idx_ups (== 1 for a trivial
  // stabilizer, the common fast path). With all-ones characters in basis_up_
  // the stored up norm is sqrt(|stabilizer|), so this is O(1) with no extra
  // storage.
  inline int64_t stab_size(int64_t idx_ups) const {
    double n = basis_up_.norm(idx_ups);
    return (int64_t)std::llround(n * n);
  }
  // All symmetries mapping `ups` to its representative (the up-stabilizer
  // coset). Unlike dns_for_ups_rep this returns an owning vector, not a span:
  // the coset of an arbitrary (generally non-representative) `ups` is computed
  // on the fly and has no persistent backing store -- pre-storing cosets for
  // every raw ups is exactly the O(raw_ups_size) memory the coupled scheme
  // avoids. It is only ever called on the non-trivial-stabilizer path (the
  // trivial path uses the single rep symmetry), so the allocation is rare.
  std::vector<int64_t> syms_ups(bit_t ups) const;

  inline gsl::span<bit_t const> dns_for_ups_rep(int64_t idx_ups) const {
    return dns_for_ups_rep_[idx_ups];
  }
  inline gsl::span<double const> norms_for_ups_rep(int64_t idx_ups) const {
    return norms_for_ups_rep_[idx_ups];
  }

  // O(1) fermi signs of a single symmetry acting on an up / dn configuration,
  // read from the precomputed tables.
  inline bool fermi_bool_ups(int64_t sym, bit_t ups) const {
    return fermi_up_.sign(sym, ups);
  }
  inline bool fermi_bool_dns(int64_t sym, bit_t dns) const {
    return fermi_dn_.sign(sym, dns);
  }

  // Block index of `dns` in up representative `idx_up`'s dn block. For electron
  // the trivial-stabilizer block is the full dn enumeration, so this is the O(1)
  // basis_dn().index() and the idx_up / dnss arguments are unused; they exist so
  // the symmetric kernels are shared verbatim with BasistJSymmetric, whose block
  // is up-rep-specific. Electron states are never absent here (no double
  // occupancy constraint), so the result is always >= 0.
  inline int64_t index_dns(bit_t dns, int64_t idx_up,
                           gsl::span<bit_t const> dnss) const {
    (void)idx_up;
    (void)dnss;
    return basis_dn_.index(dns);
  }

  // dn index and fermi sign of `dns` mapped by the single symmetry `sym`
  // (the trivial up-stabilizer path: dn block is the full enumeration). idx_up /
  // dnss are unused for electron (see index_dns).
  inline std::pair<int64_t, bool>
  index_dns_fermi(bit_t dns, int64_t sym, int64_t idx_up,
                  gsl::span<bit_t const> dnss) const {
    return {index_dns(action_.apply(sym, dns), idx_up, dnss),
            fermi_bool_dns(sym, dns)};
  }

  // dn index, fermi sign and chosen symmetry for the non-trivial up-stabilizer.
  // `syms` is the up-stabilizer coset; `dnss` is the (ascending) dn-rep block.
  // idx_dn == -1 means `dns` maps to a zero-norm (absent) dn representative.
  // Defined out-of-line: it is not a one-liner and would otherwise be inlined
  // into every kernel instantiation.
  std::tuple<int64_t, bool, int64_t>
  index_dns_fermi_sym(bit_t dns, std::vector<int64_t> const &syms,
                      gsl::span<bit_t const> dnss) const;

  // Unified "symmetrise an arbitrary (ups, dns)" primitive, the electron
  // analogue of BasisSymmetric::representative_data_fermi. Returns
  // {index, sym, norm_out, fermi_total}; index == -1 means the state has zero
  // norm (not in the basis). fermi_total already combines the up and dn fermi
  // signs of the chosen symmetry. The hot kernels do not call this per element
  // (they hoist the up-sector work); it backs index() and one-off queries, so
  // it is defined out-of-line rather than inlined into every instantiation.
  std::tuple<int64_t, int64_t, double, bool>
  representative_data_fermi(bit_t ups, bit_t dns) const;

private:
  // up sector reduced to orbit representatives (trivial all-ones characters)
  BasisSymmetric<enumeration_t> basis_up_;
  // full dn enumeration (index / iteration / size for the front block)
  BasisOnTheFly<enumeration_t> basis_dn_;

  Vector characters_;                  // the actual 1-D irrep characters
  symmetries::SitePermutation action_; // owns a copy of the group

  // O(1) fermi signs of any single symmetry acting on an up / dn configuration,
  // over the full up / dn enumerations.
  symmetries::FermiTable<enumeration_t> fermi_up_;
  symmetries::FermiTable<enumeration_t> fermi_dn_;

  // dn-rep blocks: a shared front of all dn states (norm 1) for trivial
  // up-stabilizers, followed by bespoke blocks for non-trivial ones. The
  // per-rep spans view into the flat storage and are materialised once, after
  // the storage is finalised (so they cannot dangle on reallocation).
  std::vector<bit_t> dns_storage_;
  std::vector<double> norms_storage_;
  std::vector<gsl::span<bit_t const>> dns_for_ups_rep_;
  std::vector<gsl::span<double const>> norms_for_ups_rep_;
  std::vector<int64_t> ups_offset_;

  int64_t size_ = 0;
};

// Iterates (ups_rep, dns_rep) pairs: up representatives outer, dn block inner.
template <typename enumeration_tt> class BasisElectronSymmetricIterator {
public:
  using enumeration_t = enumeration_tt;
  using bit_t = typename enumeration_t::bit_t;

  BasisElectronSymmetricIterator() = default;
  BasisElectronSymmetricIterator(
      BasisElectronSymmetric<enumeration_t> const &basis, bool begin);
  BasisElectronSymmetricIterator &operator++();
  std::pair<bit_t, bit_t> operator*() const;
  bool
  operator!=(BasisElectronSymmetricIterator<enumeration_t> const &rhs) const;

private:
  BasisElectronSymmetric<enumeration_t> const *basis_ = nullptr;
  int64_t up_idx_ = 0;
  int64_t dn_idx_ = 0;
  gsl::span<bit_t const> dnss_;
};

} // namespace xdiag::basis
