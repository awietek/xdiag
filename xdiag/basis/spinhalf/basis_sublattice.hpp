// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/extern/flat_hash_map.hpp>
#include <xdiag/extern/gsl/span>

#include <xdiag/common.hpp>

#include <xdiag/symmetries/group_action/group_action_sublattice.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/representation.hpp>

namespace xdiag::basis::spinhalf {

constexpr int64_t maximum_prefix_bits = 32;

template <typename bit_tt, int n_sublat> class BasisSublattice {
public:
  using bit_t = bit_tt;
  using iterator_t = typename std::vector<bit_t>::const_iterator;

  BasisSublattice() = default;
  BasisSublattice(Representation const &irrep);
  BasisSublattice(int64_t nup, Representation const &irrep);

  iterator_t begin() const;
  iterator_t end() const;
  int64_t size() const;
  int64_t dim() const;

  int64_t index(bit_t state) const;
  inline bit_t state(int64_t idx) const { return reps_[idx]; }
  inline double norm(int64_t idx) const { return norms_[idx]; }

  int64_t nsites() const;
  int64_t nup() const;

  Representation const &irrep() const;
  GroupActionSublattice<bit_t, n_sublat> const &group_action() const;

  bool operator==(BasisSublattice<bit_t, n_sublat> const &rhs) const;
  bool operator!=(BasisSublattice<bit_t, n_sublat> const &rhs) const;

private:
  int64_t nsites_;
  int64_t nup_;
  int64_t n_postfix_bits_;

  Representation irrep_;
  GroupActionSublattice<bit_t, n_sublat> group_action_;
  std::vector<bit_t> reps_;
  std::vector<double> norms_;
  // std::unordered_map<bit_t, gsl::span<bit_t const>> rep_search_range_;
  ska::flat_hash_map<bit_t, gsl::span<bit_t const>> rep_search_range_;

  int64_t index_of_representative(bit_t rep) const;

  // functions used in implementation of terms
public:
  inline bit_t representative(bit_t state) const {
    return group_action_.representative(state);
  }
  std::pair<int64_t, int64_t> index_sym(bit_t raw_state) const;
  std::pair<int64_t, gsl::span<int64_t const>>
  index_syms(bit_t raw_state) const;
};

} // namespace xdiag::basis::spinhalf
