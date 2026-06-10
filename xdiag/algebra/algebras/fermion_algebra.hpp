// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

#include <xdiag/algebra/algebras/algebra.hpp>

namespace xdiag::algebra {

// Spinless-fermion algebra used for symmetry analysis. Unlike
// fermion_implementation_algebra it protects no named operators: every
// compound operator (Hop, HopAsym, N, NN) is reduced to the elementary
// Cdag/C generators and the result is normal-ordered.
Algebra fermion_algebra(int64_t nsites);

} // namespace xdiag::algebra
