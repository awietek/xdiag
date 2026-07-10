// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/representation.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

XDIAG_API PermutationGroup cyclic_group(int64_t n);
XDIAG_API Representation cyclic_group_irrep(int64_t n, int64_t k);

} // namespace xdiag
