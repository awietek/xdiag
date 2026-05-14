// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/algebra/rewrite.hpp>

namespace xdiag::operators {

// MonomialRule: a size-1 SdotS{i,j} monomial expands to Exchange{i,j} +
// SzSz{i,j}. Used by spinhalf_implementation_algebra to keep those two types as
// named operators rather than collapsing them to a generic Matrix.
MonomialRule spinhalf_sdots_rule();

} // namespace xdiag::operators
