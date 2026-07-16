// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/symmetries/action/site_permutation.hpp>

namespace xdiag::symmetries {

// determines whether a state is a representative
template <typename bit_t>
bool isrepresentative(bit_t state, SitePermutation const &action);

} // namespace xdiag::symmetries
