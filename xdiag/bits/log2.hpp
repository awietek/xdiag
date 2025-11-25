// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

namespace xdiag::bits {
constexpr size_t floorlog2(size_t x) {
  return x == 1 ? 0 : 1 + floorlog2(x >> 1);
}

constexpr unsigned ceillog2(unsigned x) {
  return x == 1 ? 0 : floorlog2(x - 1) + 1;
}

} // namespace xdiag::bits
