// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

#include <xdiag/algebra/algebras/algebra.hpp>

namespace xdiag::algebra {

// Implementation algebra for the Electron block. Unlike electron_algebra (which
// reduces everything to the elementary Cdagup/Cup/Cdagdn/Cdn for the symmetry
// analysis), this keeps the named operator types that have a dedicated matrix
// kernel (see matrices/blocks/electron) and only reduces the convenience types
// to them:
//   Hop      -> Hopup + Hopdn
//   HopAsym  -> HopupAsym + HopdnAsym
//   Ntot     -> Nup + Ndn
//   Sz       -> 1/2 Nup - 1/2 Ndn
//   TotalN   -> sum_i Ntot{i}
Algebra electron_implementation_algebra(int64_t nsites);

} // namespace xdiag::algebra
