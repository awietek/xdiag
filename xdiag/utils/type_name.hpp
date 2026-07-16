// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <iostream>
#include <string_view>

namespace xdiag::utils {

// ------------------- Get compiler function string -------------------
template <typename T> constexpr std::string_view type_name_impl() {
#ifdef __clang__
  return __PRETTY_FUNCTION__;
#elif defined(__GNUC__)
  return __PRETTY_FUNCTION__;
#elif defined(_MSC_VER)
  return __FUNCSIG__;
#else
  return "unknown";
#endif
}

// ------------------- Extract the actual type name -------------------
constexpr std::string_view extract_type_name(std::string_view func_str) {
  // Find the start of the type
  std::size_t start = func_str.find("T = ");
  if (start == std::string_view::npos)
    return "unknown";
  start += 4; // skip "T = "

  // Match brackets to handle nested templates
  int depth = 0;
  std::size_t end = start;
  for (; end < func_str.size(); ++end) {
    char c = func_str[end];
    if (c == '<')
      depth++;
    if (c == '>')
      depth--;
    if (depth == 0 && (c == ']' || c == ';' || c == ' '))
      break; // stop at closing
  }

  return func_str.substr(start, end - start);
}

// ------------------- Public API -------------------
template <typename T> constexpr std::string_view get_type_name() {
  return extract_type_name(type_name_impl<T>());
}

} // namespace xdiag::utils
