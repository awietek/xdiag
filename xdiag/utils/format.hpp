// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <locale>
#include <sstream>
#include <string>

namespace xdiag::utils {

template <class T> std::string FormatWithCommas(T value) {
  std::stringstream ss;
  ss.imbue(std::locale(""));
  ss << std::fixed << value;
  return ss.str();
}

} // namespace xdiag::utils
