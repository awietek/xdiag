// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <limits>

namespace xdiag::bits {

template <typename bit_t> constexpr bit_t bitmask(int64_t length) {
  return (length == std::numeric_limits<bit_t>::digits)
             ? std::numeric_limits<bit_t>::max()
             : (((bit_t)1 << length) - 1);
}

} // namespace xdiag::bits
