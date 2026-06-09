// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

#include <xdiag/algebra/algebras/algebra.hpp>

namespace xdiag::algebra {

Algebra electron_algebra(int64_t nsites);

} // namespace xdiag::algebra
