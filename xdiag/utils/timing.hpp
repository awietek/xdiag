// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <chrono>
#include <string>

#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

using clock_t = std::chrono::high_resolution_clock;

XDIAG_API auto rightnow() -> decltype(clock_t::now());
XDIAG_API void timing(std::chrono::time_point<clock_t> const &t0,
                      std::chrono::time_point<clock_t> const &t1,
                      std::string msg = "", int verbosity = 0);
XDIAG_API void tic(bool begin = true, std::string msg = "", int verbosity = 0);
XDIAG_API void toc(std::string msg = "", int verbosity = 0);

} // namespace xdiag
