// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <string_view>
#include <vector>

#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/utils/type_name.hpp>

namespace xdiag::symmetries {

// Precomputed table of fermi signs over a state enumeration and a permutation
// group: sign(sym, state) == fermi_bool_of_permutation(state, group[sym]).
//
// fermi_bool_of_permutation is O(nsites); evaluating it per matrix element in
// the symmetric apply/matrix kernels would dominate the hot loop. This trades
// |group| x |enumeration| bits of one-time storage for an O(1) lookup (plus the
// enumeration's own index(), O(1) for the LinTable enumerations the symmetric
// blocks use). Templated on the enumeration, like RepresentativeTable.
template <typename enumeration_tt> class FermiTable {
public:
  using enumeration_t = enumeration_tt;
  using bit_t = typename enumeration_t::bit_t;
  static constexpr std::string_view type_name =
      utils::get_type_name<FermiTable<enumeration_t>>();

  FermiTable() = default;
  FermiTable(enumeration_t const &enumeration, PermutationGroup const &group);

  inline bool sign(int64_t sym, bit_t state) const {
    return table_[sym * size_ + enumeration_.index(state)];
  }

private:
  enumeration_t enumeration_;
  int64_t size_ = 0; // number of states in the enumeration
  std::vector<bool> table_;
};

} // namespace xdiag::symmetries
