// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <string_view>

#include <xdiag/basis/basis.hpp>
#include <xdiag/states/product_state.hpp>
#include <xdiag/utils/type_name.hpp>

namespace xdiag::basis {

template <typename enumeration_tt>
class BasisOnTheFly : public BasisType<BasisOnTheFly<enumeration_tt>> {
public:
  using enumeration_t = enumeration_tt;
  using bit_t = typename enumeration_t::bit_t;
  using iterator_t = typename enumeration_t::iterator_t;
  static constexpr std::string_view type_name =
      utils::get_type_name<BasisOnTheFly<enumeration_t>>();

  BasisOnTheFly() = default;
  BasisOnTheFly(enumeration_t const &enumeration);

  int64_t size() const override;
  int64_t nsites() const override;
  int64_t d() const; // Local Hilbert space dimension per site
  int64_t index(bit_t bits) const;
  iterator_t begin() const;
  iterator_t end() const;

  bool operator==(BasisOnTheFly<enumeration_t> const &rhs) const;
  bool operator!=(BasisOnTheFly<enumeration_t> const &rhs) const;

private:
  enumeration_t enumeration_;
};

} // namespace xdiag::basis
