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
#include <xdiag/states/product_state.hpp>
#include <xdiag/symmetries/representation.hpp>
#include <xdiag/symmetries/tables/representative_table.hpp>
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
  BasisSymmetric(Representation const &irrep, enumeration_t const &enumeration);

  int64_t size() const;
  int64_t nsites() const;
  Representation const& irrep() const;
  
  inline std::tuple<int64_t, int64_t, double>
  representative_index_symmetry_norm(bit_t bits) const {
    int64_t idx = enumeration_.index(bits);
    return {table_.representative_index(idx),
            table_.representative_symmetry(idx),
            table_.representative_norm(idx)};
  }
  iterator_t begin() const;
  iterator_t end() const;

  ProductState
  product_state(int64_t idx,
                std::vector<std::string> const &dict) const override;

  bool operator==(BasisSymmetric<enumeration_t> const &rhs) const;
  bool operator!=(BasisSymmetric<enumeration_t> const &rhs) const;

private:
  Representation irrep_;
  enumeration_t enumeration_;
  table_t table_;
};

} // namespace xdiag::basis
