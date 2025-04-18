// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <utility>
#include <vector>
#include <xdiag/extern/gsl/span>

#include <xdiag/combinatorics/subsets_indexing.hpp>
#include <xdiag/common.hpp>
#include <xdiag/symmetries/group_action/group_action_lookup.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/representation.hpp>

namespace xdiag::basis::spinhalf {

template <typename bit_tt> class BasisSymmetricNoSz {
public:
  using bit_t = bit_tt;
  using iterator_t = typename std::vector<bit_t>::const_iterator;
  using span_size_t = gsl::span<int64_t const>::size_type;

  BasisSymmetricNoSz() = default;
  BasisSymmetricNoSz(Representation const &irrep);

  int64_t dim() const;
  int64_t size() const;
  iterator_t begin() const;
  iterator_t end() const;

  inline int64_t index(bit_t state) const {
    return index_for_rep_[subsets_basis_.index(state)];
  }
  inline bit_t state(int64_t idx) const { return reps_[idx]; }
  inline double norm(int64_t idx) const { return norms_[idx]; }

  int64_t nsites() const;
  GroupActionLookup<bit_t> const &group_action() const;
  Representation const &irrep() const;

  bool operator==(BasisSymmetricNoSz const &rhs) const;
  bool operator!=(BasisSymmetricNoSz const &rhs) const;

private:
  int64_t nsites_;
  GroupActionLookup<bit_t> group_action_;
  Representation irrep_;
  combinatorics::SubsetsIndexing<bit_t> subsets_basis_;

  std::vector<bit_t> reps_;
  std::vector<int64_t> index_for_rep_;
  std::vector<int64_t> syms_;
  std::vector<std::pair<span_size_t, span_size_t>> sym_limits_for_rep_;
  std::vector<double> norms_;

  int64_t size_;

  // functions used in implementation of terms
public:
  inline bit_t representative(bit_t raw_state) const {
    return reps_[index(raw_state)];
  }

  inline std::pair<int64_t, int64_t> index_sym(bit_t raw_state) const {
    int64_t raw_idx = subsets_basis_.index(raw_state);
    int64_t index = index_for_rep_[raw_idx];
    if (index == invalid_index) {
      return {invalid_index, 0};
    }
    int64_t start = sym_limits_for_rep_[raw_idx].first;
    return {index, syms_[start]};
  }

  inline std::pair<int64_t, gsl::span<int64_t const>>
  index_syms(bit_t raw_state) const {
    int64_t raw_idx = subsets_basis_.index(raw_state);
    int64_t index = index_for_rep_[raw_idx];
    auto [start, length] = sym_limits_for_rep_[raw_idx];
    return {index, {syms_.data() + start, length}};
  }
};

} // namespace xdiag::basis::spinhalf
