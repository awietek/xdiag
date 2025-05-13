// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/common.hpp>

#include <string>

namespace xdiag {

XDIAG_API std::string version_string();
XDIAG_API void say_hello();
XDIAG_API void print_version();

} // namespace xdiag
