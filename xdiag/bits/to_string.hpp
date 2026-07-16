// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <cstdint>
#include <limits>
#include <string>

namespace xdiag::bits {

template <typename bit_t>
std::string to_string(bit_t bits,
                      int64_t size = std::numeric_limits<bit_t>::digits,
                      bool reverse = true);
}
