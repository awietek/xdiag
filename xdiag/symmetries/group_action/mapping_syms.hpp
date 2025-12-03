// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <vector>

namespace xdiag::symmetries {

// Computes the symmetries, which map origin to target
template <typename bit_t, class group_action_t>
std::vector<int64_t> mapping_syms(bit_t origin, bit_t target,
                                  group_action_t const &group_action);

} // namespace xdiag::symmetries
