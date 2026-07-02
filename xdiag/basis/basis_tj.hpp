// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include <xdiag/basis/basis.hpp>
#include <xdiag/basis/basis_onthefly.hpp>
#include <xdiag/bits/bitmask.hpp>
#include <xdiag/bits/extract_deposit.hpp>
#include <xdiag/bits/get_set.hpp>
#include <xdiag/bits/popcount.hpp>
#include <xdiag/bits/zero_one.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/utils/type_name.hpp>

namespace xdiag::basis {

template <typename enumeration_tt> class BasistJIterator;

// Spinful tJ basis (local dimension d = 3: empty, up, dn) without permutation
// symmetry. Like BasisElectron it splits into an "ups" sector and a "dns"
// sector, but the no-double-occupancy constraint (ups & dns == 0) couples them:
// a dn can only sit on a site not occupied by an up. We store the dn sector
// COMPRESSED into the nsites - nup complementary sites:
//
//   ups   in enumeration on (nsites, nup)         -- the up configuration
//   dncs  = the dn config packed (in site order) into the non-up sites
//
// The linear index is  ups_offset(idx_up) + index_dncs(ups, dncs).
//
// Number conservation vs not is exactly the choice of enumeration (as in the
// electron block), and it is confined to a small "seam" the kernels go through
// -- ups_offset(idx_up), the compressed dn enumeration for a given ups, and
// index_dncs(ups, dncs):
//
//   np     : ups in Combinations/LinTable(nsites,nup); the dn space is
//            LinTable(nsites-nup, ndn), the SAME for every ups, so
//            ups_offset(idx_up) = idx_up * size_dncs and index_dncs is a
//            LinTable lookup.
//   no-np  : ups in Subsets(nsites); the dn space is Subsets(nsites-popcount(ups)),
//            which VARIES per ups. Because that fiber depends only on
//            popcount(ups), there are only nsites+1 distinct fibers; we store
//            them (cheap Subsets) and index by popcount(ups), so basis_dncs(ups)
//            is still a const-ref (no per-ups table copy). ups_offset is a stored
//            prefix sum and index_dncs is the value itself.
//
// The option-B kernels operate in this compressed dn space directly -- a full
// site s maps to the compressed rank popcount((~ups) & below(s)), computed once
// per ups -- so no pdep/pext is needed on the hot path, and the same kernels
// serve np and no-np through the seam below. The full dn configuration is only
// reconstructed when materialising a ProductState (iteration) or on the cold
// Cdag/C-string fallback, via tj_decompress_dns -- never by the named kernels.
template <typename enumeration_tt>
class BasistJ : public BasisType<BasistJ<enumeration_tt>> {
public:
  using enumeration_t = enumeration_tt;
  using bit_t = typename enumeration_t::bit_t;
  using iterator_t = BasistJIterator<enumeration_t>;
  static constexpr std::string_view type_name =
      utils::get_type_name<BasistJ<enumeration_t>>();

  BasistJ() = default;
  // np: ups on (nsites,nup), compressed dn on (nsites-nup,ndn).
  BasistJ(enumeration_t const &enum_up, enumeration_t const &enum_dncs);
  // no-np: ups over all subsets of nsites; the dn fibers are derived (Subsets on
  // the per-ups complement). Only valid for enumeration_t == Subsets.
  explicit BasistJ(enumeration_t const &enum_up);

  int64_t size() const override;
  int64_t nsites() const override;
  int64_t index(ProductState const &pstate) const override;
  constexpr int64_t d() const { return 3; }

  iterator_t begin() const;
  iterator_t end() const;

  // ---- kernel-facing seam (np/no-np-agnostic) --------------------------------
  // The up sector (full nsites) -- iteration, index_up and size, reached through
  // directly by the kernels (as in the electron block).
  BasisOnTheFly<enumeration_t> const &basis_up() const { return basis_up_; }

  // The compressed dn enumeration to iterate for a given ups. np stores a single
  // shared fiber; no-np stores one fiber per popcount(ups). Either way this is a
  // const-ref to a stored fiber, so the hot loop never copies a table.
  inline BasisOnTheFly<enumeration_t> const &basis_dncs(bit_t ups) const {
    return (dncs_.size() == 1) ? dncs_[0] : dncs_[bits::popcount(ups)];
  }
  // Linear offset of the dn block of up configuration idx_up: a stored prefix sum
  // so the same accessor serves np (linear) and no-np (varying fiber).
  inline int64_t ups_offset(int64_t idx_up) const {
    return ups_offset_[idx_up];
  }
  // Index of a compressed dn configuration within its ups block.
  inline int64_t index_dncs(bit_t ups, bit_t dncs) const {
    return basis_dncs(ups).index(dncs);
  }
  // Index of an up configuration in the up sector (needed by up-changing
  // kernels to look up the offset of the target ups block).
  inline int64_t index_up(bit_t ups) const { return basis_up_.index(ups); }

  bool operator==(BasistJ<enumeration_t> const &rhs) const;
  bool operator!=(BasistJ<enumeration_t> const &rhs) const;

private:
  BasisOnTheFly<enumeration_t> basis_up_; // ups (np: (nsites,nup); no-np: subsets)
  // Compressed dn fibers: size 1 for np (the shared fiber), size nsites+1 for
  // no-np (indexed by popcount(ups), the complement Subsets(nsites-popcount)).
  std::vector<BasisOnTheFly<enumeration_t>> dncs_;
  std::vector<int64_t> ups_offset_; // prefix sum of per-ups dn sizes
  int64_t size_ = 0;
};

// Decompress a compressed dn configuration `dncs` (on the nsites-nup non-up
// sites) back to a full nsites dn configuration, given `ups`: a parallel bits
// deposit of `dncs` into the non-up sites (pdep on BMI2, the Bitset overload of
// bits::deposit for BitsetDynamic). Used off the option-B hot path -- by
// ProductState iteration and the general Cdag/C-string fallback, never by the
// named kernels (which keep the dn compressed and never materialise it).
template <typename bit_t>
inline bit_t tj_decompress_dns(bit_t ups, bit_t dncs, int64_t nsites) {
  bit_t complement = (~ups) & bits::bitmask<bit_t>(nsites, nsites);
  return bits::deposit(dncs, complement);
}

// Compress a full nsites dn configuration `dns` into the nsites-nup non-up sites
// of `ups` (the inverse of tj_decompress_dns): a parallel bits extract (pext on
// BMI2, the Bitset overload of bits::extract for BitsetDynamic). `dns` is assumed
// to have no bit on an up site (the no-double-occupancy invariant).
template <typename bit_t>
inline bit_t tj_compress_dns(bit_t ups, bit_t dns, int64_t nsites) {
  bit_t complement = (~ups) & bits::bitmask<bit_t>(nsites, nsites);
  return bits::extract(dns, complement);
}

// Iterates (ups, dns) pairs with dns DECOMPRESSED to full nsites, so the generic
// local_state(pair, i) works: up sector outer, compressed dn sector inner. The
// dn fiber is re-fetched from the basis when ups advances, so the same iterator
// serves np (constant fiber) and no-np (fiber depends on popcount(ups)).
template <typename enumeration_tt> class BasistJIterator {
public:
  using enumeration_t = enumeration_tt;
  using bit_t = typename enumeration_t::bit_t;
  using sub_iterator_t = typename enumeration_t::iterator_t;

  BasistJIterator() = default;
  BasistJIterator(BasistJ<enumeration_t> const &basis, bool begin);
  BasistJIterator &operator++();
  std::pair<bit_t, bit_t> operator*() const; // {ups, full dns}
  bool operator!=(BasistJIterator<enumeration_t> const &rhs) const;

private:
  BasistJ<enumeration_t> const *basis_ = nullptr;
  sub_iterator_t up_, up_end_;
  sub_iterator_t dn_, dn_end_;
  int64_t nsites_ = 0;
};

} // namespace xdiag::basis
