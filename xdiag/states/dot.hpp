// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/math/complex.hpp>
#include <xdiag/states/state.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

XDIAG_API double dot(State const &v, State const &w);
XDIAG_API complex dotC(State const &v, State const &w);

XDIAG_API arma::mat matrix_dot(State const &v, State const &w);
XDIAG_API arma::cx_mat matrix_dotC(State const &v, State const &w);

} // namespace xdiag
