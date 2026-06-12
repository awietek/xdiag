// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

#include <xdiag/basis/basis.hpp>
#include <xdiag/math/vector.hpp>
#include <xdiag/states/product_state.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/tables/representative_table.hpp>
#include <xdiag/utils/likely.hpp>
#include <xdiag/utils/type_name.hpp>

namespace xdiag::basis {

template <typename enumeration_tt>
class BasisSymmetric : public BasisType<BasisSymmetric<enumeration_tt>> {
public:
  using enumeration_t = enumeration_tt;
  using bit_t = typename enumeration_t::bit_t;
  using table_t = symmetries::RepresentativeTable<enumeration_t>;
  using iterator_t = typename table_t::const_iterator;
  static constexpr std::string_view type_name =
      utils::get_type_name<BasisSymmetric<enumeration_t>>();

  BasisSymmetric() = default;
  BasisSymmetric(enumeration_t const &enumeration,
                 PermutationGroup const &group, Vector const &characters,
                 bool fermionic = false);

  int64_t size() const override;
  int64_t nsites() const;
  int64_t d() const; // Local Hilbert space dimension per site
  PermutationGroup const &group() const;
  Vector const &characters() const;

  // Returns {raw_rep_idx, sym, norm_out}.
  // raw_rep_idx == 0 means zero-norm (invalid); actual index is raw_rep_idx
  // - 1.
  inline std::tuple<int64_t, int64_t, double>
  representative_data(bit_t bits) const {
    int64_t idx = enumeration_.index(bits);
    int64_t raw = table_.raw_representative_index(idx);
    if (XDIAG_LIKELY(raw)) {
      return {raw, table_.representative_symmetry(idx),
              table_.representative_norm(raw - 1)};
    } else {
      return {0, 0, 0.0};
    }
  }
  // Like representative_data, but also returns whether the mapping to the
  // representative picks up a fermi minus sign (true => -1). Used by
  // term_offdiag_fermionic; the bosonic term_offdiag_sym keeps using the
  // 3-tuple representative_data.
  inline std::tuple<int64_t, int64_t, double, bool>
  representative_data_fermi(bit_t bits) const {
    int64_t idx = enumeration_.index(bits);
    int64_t raw = table_.raw_representative_index(idx);
    if (XDIAG_LIKELY(raw)) {
      return {raw, table_.representative_symmetry(idx),
              table_.representative_norm(raw - 1),
              table_.representative_fermi_bool(idx)};
    } else {
      return {0, 0, 0.0, false};
    }
  }
  inline double norm(int64_t idx) const {
    return table_.representative_norm(idx);
  }
  inline double inv_norm(int64_t idx) const {
    return table_.inv_representative_norm(idx);
  }

  inline bit_t operator[](int64_t idx) const { return table_[idx]; }
  iterator_t begin() const;
  iterator_t end() const;

  bool operator==(BasisSymmetric<enumeration_t> const &rhs) const;
  bool operator!=(BasisSymmetric<enumeration_t> const &rhs) const;

private:
  enumeration_t enumeration_;
  PermutationGroup group_;
  Vector characters_;
  table_t table_;
};

} // namespace xdiag::basis
