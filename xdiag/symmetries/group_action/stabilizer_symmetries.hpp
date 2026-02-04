// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <vector>

namespace xdiag::symmetries {

// Computes stabilizing symmetries of "bits" when group is applied
template <typename bit_t, class group_action_t>
std::vector<int64_t> stabilizer_symmetries(bit_t bits,
                                           group_action_t const &group_action);
} // namespace xdiag::symmetries
