// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <sstream>
#include <string>

namespace xdiag::utils {

template <typename T> inline std::string to_string_generic(T const &x) {
  std::ostringstream ss;
  ss << x;
  return ss.str();
}

} // namespace xdiag::utils
