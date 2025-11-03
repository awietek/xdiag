// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "basis_tj_distributed.hpp"

namespace xdiag::basis {

int64_t dim(BasistJDistributed const &basis) {
  return std::visit([&](auto &&b) { return b.dim(); }, basis);
}
int64_t size(BasistJDistributed const &basis) {
  return std::visit([&](auto &&b) { return b.size(); }, basis);
}
int64_t size_max(BasistJDistributed const &basis) {
  return std::visit([&](auto &&b) { return b.size_max(); }, basis);
}
int64_t size_min(BasistJDistributed const &basis) {
  return std::visit([&](auto &&b) { return b.size_min(); }, basis);
}

template <typename bit_t> bool has_bit_t(BasistJDistributed const &basis) try {
  return std::visit(
      [](auto &&b) {
        using basis_t = typename std::decay<decltype(b)>::type;
        return std::is_same<bit_t, typename basis_t::bit_t>::value;
      },
      basis);
}
XDIAG_CATCH

template bool has_bit_t<uint32_t>(BasistJDistributed const &basis);
template bool has_bit_t<uint64_t>(BasistJDistributed const &basis);

} // namespace xdiag::basis
