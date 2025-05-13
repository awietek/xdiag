// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "basis_electron_distributed.hpp"

namespace xdiag::basis {

int64_t dim(BasisElectronDistributed const &basis) {
  return std::visit([&](auto &&b) { return b.dim(); }, basis);
}
int64_t size(BasisElectronDistributed const &basis) {
  return std::visit([&](auto &&b) { return b.size(); }, basis);
}
int64_t size_max(BasisElectronDistributed const &basis) {
  return std::visit([&](auto &&b) { return b.size_max(); }, basis);
}
int64_t size_min(BasisElectronDistributed const &basis) {
  return std::visit([&](auto &&b) { return b.size_min(); }, basis);
}

template <typename bit_t>
bool has_bit_t(BasisElectronDistributed const &basis) try {
  return std::visit(
      [](auto &&b) {
        using basis_t = typename std::decay<decltype(b)>::type;
        return std::is_same<bit_t, typename basis_t::bit_t>::value;
      },
      basis);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}
template bool has_bit_t<uint32_t>(BasisElectronDistributed const &basis);
template bool has_bit_t<uint64_t>(BasisElectronDistributed const &basis);

} // namespace xdiag::basis
