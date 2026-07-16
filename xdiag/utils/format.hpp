// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <complex>
#include <locale>

// fmt in compiled mode: core defined once in utils/format.cpp, so
// FMT_HEADER_ONLY must not be defined anywhere (inline vs compiled defs clash).
#include <extern/fmt/format.hpp>

namespace fmt {

// Formatter for complex numbers
template <> struct formatter<std::complex<double>> {
  template <typename ParseContext> constexpr auto parse(ParseContext &ctx) {
    return ctx.begin();
  }

  template <typename FormatContext>
  inline auto format(std::complex<double> const &number, FormatContext &ctx) {
    if (std::imag(number) < 0.) {
      return fmt::format_to(ctx.out(), "{0}-i{1}", std::real(number),
                            -std::imag(number));
    } else {
      return fmt::format_to(ctx.out(), "{0}+i{1}", std::real(number),
                            std::imag(number));
    }
  }
};

// comma separator for large numbers after 1,000 and 1,000,000 etc
struct german_punct : std::numpunct<char> {
protected:
  char do_decimal_point() const override { return ','; }
  char do_thousands_sep() const override { return '.'; }
  std::string do_grouping() const override { return "\3"; }
};

template <typename... Args>
std::string format_de(fmt::format_string<Args...> fmt_str, Args &&...args) {
  static std::locale german(std::locale::classic(), new german_punct);

  return fmt::vformat(
      german,
      fmt::string_view(fmt_str), // ← THIS is the important trick
      fmt::make_format_args(args...));
}

} // namespace fmt
