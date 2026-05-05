// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <utility>
#include <xdiag/extern/gsl/span>
#include <xdiag/symmetries/action/site_permutation.hpp>

namespace xdiag::symmetries {

template <typename bit_t>
bit_t representative(bit_t state, SitePermutation const &action);

template <typename bit_t>
std::pair<bit_t, int64_t> representative_sym(bit_t state,
                                             SitePermutation const &action);

template <typename bit_t>
std::pair<bit_t, std::vector<int64_t>>
representative_syms(bit_t state, SitePermutation const &action);

} // namespace xdiag::symmetries
