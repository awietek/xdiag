// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/armadillo.hpp>
#include <xdiag/states/state.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

// Various norms
XDIAG_API double norm(State const &v);
XDIAG_API double norm1(State const &v);
XDIAG_API double norminf(State const &v);

} // namespace xdiag
