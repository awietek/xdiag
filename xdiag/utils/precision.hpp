#pragma once

#include <xdiag/common.hpp>

namespace xdiag::numeric {

// Predifined absolute precision for comparisions
static constexpr float float_atol =
    10000. * std::numeric_limits<float>::epsilon();
static constexpr double double_atol =
    10000. * std::numeric_limits<double>::epsilon();

template <class T> struct atol {
  static constexpr real_t<T> val();
};

template <> struct atol<float> {
  static constexpr real_t<float> val() { return float_atol; }
};

template <> struct atol<double> {
  static constexpr real_t<double> val() { return double_atol; }
};

template <> struct atol<std::complex<float>> {
  static constexpr real_t<std::complex<float>> val() { return float_atol; }
};

template <> struct atol<std::complex<double>> {
  static constexpr real_t<std::complex<double>> val() { return double_atol; }
};

// Predifined relative precision for comparisions
static constexpr float float_rtol =
    100000. * std::numeric_limits<float>::epsilon();
static constexpr double double_rtol =
    100000. * std::numeric_limits<double>::epsilon();

template <class T> struct rtol {
  static constexpr real_t<T> val();
};

template <> struct rtol<float> {
  static constexpr real_t<float> val() { return float_rtol; }
};

template <> struct rtol<double> {
  static constexpr real_t<double> val() { return double_rtol; }
};

template <> struct rtol<std::complex<float>> {
  static constexpr real_t<std::complex<float>> val() { return float_rtol; }
};

template <> struct rtol<std::complex<double>> {
  static constexpr real_t<std::complex<double>> val() { return double_rtol; }
};

  
} // namespace xdiag
