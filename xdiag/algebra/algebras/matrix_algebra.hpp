// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

#include <xdiag/algebra/algebras/algebra.hpp>

namespace xdiag::algebra {

// `d` is the local Hilbert space dimension per site (matrices passed in via
// "Matrix" Ops must be d^nsites in each dimension).
Algebra matrix_algebra(int64_t nsites, int64_t d);

} // namespace xdiag::algebra
