// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <limits>

namespace xdiag::bits {

template <typename bit_t, typename length_t>
inline bit_t bitmask(length_t length) {
  return (length == std::numeric_limits<bit_t>::digits)
             ? std::numeric_limits<bit_t>::max()
             : (((bit_t)1 << length) - 1);
}

} // namespace xdiag::bits
