// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

#include <xdiag/algebra/rewrite.hpp>

namespace xdiag::operators {

// MonomialRule: combine the first pair of adjacent "Matrix" ops (both with
// sites) into a single "Matrix" op using combine_matrix_ops. `d` is the local
// Hilbert space dimension per site.
MonomialRule combine_matrix_rule(int64_t d);

} // namespace xdiag::operators
