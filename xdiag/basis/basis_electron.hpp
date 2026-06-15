// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <memory>
#include <string_view>
#include <utility>

#include <xdiag/basis/basis.hpp>
#include <xdiag/basis/basis_onthefly.hpp>
#include <xdiag/states/product_state.hpp>
#include <xdiag/utils/product_iterator.hpp>
#include <xdiag/utils/type_name.hpp>

namespace xdiag::basis {

// Spinful electron basis (local dimension d = 4: empty, up, dn, up&dn) without
// permutation symmetry. It is the tensor product of an "ups" and a "dns"
// enumeration of the same type (e.g. both LinTable<uint32_t> with different
// particle numbers, or both Subsets<uint32_t> for no number conservation).
// A basis state is the pair (ups, dns); the linear index is
// index_up(ups) * size_dn + index_dn(dns) (up sector outer, dn sector inner).
template <typename enumeration_tt>
class BasisElectron : public BasisType<BasisElectron<enumeration_tt>> {
public:
  using enumeration_t = enumeration_tt;
  using bit_t = typename enumeration_t::bit_t;
  // Iterates (ups, dns) pairs: up sector outer, dn sector inner.
  using iterator_t =
      utils::product_iterator<typename enumeration_t::iterator_t,
                              typename enumeration_t::iterator_t>;
  static constexpr std::string_view type_name =
      utils::get_type_name<BasisElectron<enumeration_t>>();

  BasisElectron() = default;
  BasisElectron(enumeration_t const &enum_up, enumeration_t const &enum_dn);

  int64_t size() const override;   // size_up * size_dn
  int64_t nsites() const override; // sites of the up (== dn) enumeration
  constexpr int64_t d() const { return 4; }

  // Linear index of the basis state (ups, dns).
  int64_t index(bit_t ups, bit_t dns) const;

  iterator_t begin() const;
  iterator_t end() const;

  BasisOnTheFly<enumeration_t> const &basis_up() const;
  BasisOnTheFly<enumeration_t> const &basis_dn() const;

  bool operator==(BasisElectron<enumeration_t> const &rhs) const;
  bool operator!=(BasisElectron<enumeration_t> const &rhs) const;

private:
  BasisOnTheFly<enumeration_t> basis_up_;
  BasisOnTheFly<enumeration_t> basis_dn_;
};

} // namespace xdiag::basis
