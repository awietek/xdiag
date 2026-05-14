// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

#include <xdiag/algebra/rewrite.hpp>

namespace xdiag::operators {

// OpRule: if a "Matrix" op has unsorted sites, sort them and permute the
// matrix accordingly using permute_matrix_op. `d` is the local Hilbert space
// dimension per site.
OpRule sort_matrix_sites_rule(int64_t d);

} // namespace xdiag::operators
