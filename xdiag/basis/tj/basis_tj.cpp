// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "basis_tj.hpp"

namespace xdiag::basis {

int64_t dim(BasistJ const &basis) {
  return std::visit([&](auto &&b) { return b.dim(); }, basis);
}
int64_t size(BasistJ const &basis) {
  return std::visit([&](auto &&b) { return b.size(); }, basis);
}

template <typename bit_t> bool has_bit_t(BasistJ const &basis) try {
  return std::visit(
      [](auto &&b) {
        using basis_t = typename std::decay<decltype(b)>::type;
        return std::is_same<bit_t, typename basis_t::bit_t>::value;
      },
      basis);
}
XDIAG_CATCH

template bool has_bit_t<uint32_t>(BasistJ const &basis);
template bool has_bit_t<uint64_t>(BasistJ const &basis);

} // namespace xdiag::basis
