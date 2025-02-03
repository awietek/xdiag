#pragma once

#include <cassert>
#include <complex>
#include <cstdint>
#include <limits>
#include <type_traits>

#define FMT_HEADER_ONLY
#include <xdiag/extern/fmt/format.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/xdiag_api.hpp>

#ifdef XDIAG_JULIA_WRAPPER
#define XDIAG_OFFSET 1
#else
#define XDIAG_OFFSET 0
#endif

namespace xdiag {

// helper functions for complex numbers
using complex = std::complex<double>;
template <typename T> constexpr bool isreal() {
  return std::is_same<double, T>::value;
}
constexpr double real(double x) { return x; }
constexpr double real(complex x) { return x.real(); }

constexpr double imag(double) { return 0.; }
constexpr double imag(complex x) { return x.imag(); }

constexpr double conj(double x) { return x; }
inline complex conj(complex x) { return std::conj(x); }

// Comment: check if replacable by std::optional
constexpr int64_t invalid_index = (int64_t)-1;
constexpr int64_t undefined = std::numeric_limits<int64_t>::min();

// Helper type for visitor patterns
template <class... Ts> struct overload : Ts... { using Ts::operator()...; };
template <class... Ts> overload(Ts...) -> overload<Ts...>;

template <typename T> inline std::string to_string_generic(T const &x) {
  std::ostringstream ss;
  ss << x;
  return ss.str();
}

} // namespace xdiag

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

} // namespace fmt
