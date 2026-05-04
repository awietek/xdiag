// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <vector>
#include <xdiag/bits/bitvector.hpp>
#include <xdiag/symmetries/action/site_permutation.hpp>
#include <xdiag/symmetries/representation.hpp>

namespace xdiag::symmetries {

// RepresentativeTable maps every state in an enumeration to the
// lexicographically smallest state in its symmetry orbit (the "representative")
// under a given SitePermutation action and 1-D Representation (irrep).
//
// For each state at enumeration index idx, the table stores:
//   representative(idx)         – the orbit representative (a bit state)
//   index_of_representative(idx) – index of the representative in the rep list
//   symmetry(idx)               – index s such that action.apply(s, state) ==
//                                  representative(idx)
//   norm(idx)                   – sqrt(|orbit|) weight for symmetry-adapted
//                                  basis states (zero-norm states are excluded)
//
// Iteration (begin/end) yields the representatives in enumeration order.
// The group of the SitePermutation must equal the group of the Representation.
template <typename enumeration_tt> class RepresentativeTable {
public:
  using enumeration_t = enumeration_tt;
  using bit_t = typename enumeration_t::bit_t;
  using const_iterator = bits::BitVectorConstIterator<bit_t>;

  RepresentativeTable() = default;
  RepresentativeTable(Representation const &irrep,
                      enumeration_t const &enumeration);

  inline bit_t operator[](int64_t idx) const { return representative_[idx]; }
  inline bit_t representative(int64_t idx) const {
    return representative_[representative_index_[idx] - 1];
  }
  // Returns 0 if state has zero norm (invalid), actual rep index + 1 otherwise
  inline int64_t raw_representative_index(int64_t idx) const {
    return (int64_t)representative_index_[idx];
  }
  inline int64_t representative_index(int64_t idx) const {
    return (int64_t)representative_index_[idx] - 1;
  }
  inline int64_t representative_symmetry(int64_t idx) const {
    return (int64_t)representative_symmetry_[idx];
  }
  inline double representative_norm(int64_t idx) const {
    return norm_[representative_norm_index_[idx]];
  }
  int64_t size() const;

  const_iterator begin() const noexcept { return representative_.cbegin(); }
  const_iterator end() const noexcept { return representative_.cend(); }

  bool operator==(RepresentativeTable const &rhs) const;
  bool operator!=(RepresentativeTable const &rhs) const;

private:
  bits::BitVector<bit_t> representative_;
  bits::BitVector<uint64_t> representative_index_;
  bits::BitVector<uint64_t> representative_symmetry_;
  bits::BitVector<uint64_t> representative_norm_index_;
  std::vector<double> norm_;
};

} // namespace xdiag::symmetries
