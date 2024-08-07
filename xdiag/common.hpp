#pragma once

#include <cassert>
#include <complex>
#include <cstdint>
#include <math.h>
#include <utility>

#include <xdiag/bits/bitops.hpp>
#include <xdiag/config.hpp>
#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/xdiag_show.hpp>

// TODO: IS THIS REALLY NECESASARY
#define BSTR(x) bits::bits_to_string(x, n_sites)

namespace xdiag::detail {

template <class T> struct is_complex_t : public std::false_type {};
template <class T>
struct is_complex_t<std::complex<T>> : public std::true_type {};

// Complex real/imag/conj
template <class coeff_t> struct real_type_struct { typedef coeff_t type; };

template <class coeff_t> struct real_type_struct<std::complex<coeff_t>> {
  typedef coeff_t type;
};

template <class coeff_t> struct complex_type_struct {
  typedef std::complex<coeff_t> type;
};

template <class coeff_t> struct complex_type_struct<std::complex<coeff_t>> {
  typedef std::complex<coeff_t> type;
};

} // namespace xdiag::detail

namespace xdiag {

using std_bit_t = uint64_t;
using complex = std::complex<double>;

using namespace std::literals::complex_literals;
constexpr double pi = M_PI;
inline complex operator*(complex a, int b) { return a * (double)b; }
inline complex operator*(int a, complex b) { return b * a; }
inline complex operator/(complex a, int b) { return a / (double)b; }
inline complex operator/(int a, complex b) { return (double)a / b; }

template <class T> constexpr bool iscomplex() {
  return detail::is_complex_t<T>::value;
}
template <class T> constexpr bool isreal() {
  return !detail::is_complex_t<T>::value;
}

constexpr int64_t invalid_index = (int64_t)-1;
constexpr int64_t undefined = std::numeric_limits<int64_t>::min();
constexpr double default_double = std::numeric_limits<double>::max();

constexpr bool index_not_found(int64_t idx) { return idx < 0; }
constexpr bool index_valid(int64_t idx) { return idx >= 0; }

template <class coeff_t>
using real_t = typename detail::real_type_struct<coeff_t>::type;

template <class coeff_t>
using complex_t = typename detail::complex_type_struct<coeff_t>::type;

inline float real(float x) { return x; }
inline double real(double x) { return x; }
inline float real(std::complex<float> x) { return x.real(); }
inline double real(std::complex<double> x) { return x.real(); }

inline float imag(float) { return 0.; }
inline double imag(double) { return 0.; }
inline float imag(std::complex<float> x) { return x.imag(); }
inline double imag(std::complex<double> x) { return x.imag(); }

inline float conj(float x) { return x; }
inline double conj(double x) { return x; }
inline std::complex<float> conj(std::complex<float> x) { return std::conj(x); }
inline std::complex<double> conj(std::complex<double> x) {
  return std::conj(x);
}

inline arma::cx_mat to_cx_mat(arma::mat const &A) {
  return arma::cx_mat(A, arma::mat(A.n_rows, A.n_cols, arma::fill::zeros));
}

inline arma::cx_mat to_cx_mat(arma::cx_mat const &A) { return A; }

inline bool file_exists(std::string filename) {
  std::ifstream infile(filename.c_str());
  return infile.good();
}

template <typename coeff_t>
void resize_vector(std::vector<coeff_t> &vec, int64_t size) {
  try {
    vec.resize(size);
  } catch (...) {
    XDIAG_THROW("Unable to resize vector");
  }
}

// Helper type for visitor patterns
template <class... Ts> struct overload : Ts... { using Ts::operator()...; };
template <class... Ts> overload(Ts...) -> overload<Ts...>;

template <typename T> inline std::string to_string_generic(T const &x) {
  std::ostringstream ss;
  ss << x;
  return ss.str();
}

} // namespace xdiag

// Formatter for complex numbers
namespace fmt {

template <> struct formatter<std::complex<double>> {
  template <typename ParseContext> constexpr auto parse(ParseContext &ctx) {
    return ctx.begin();
  }

  template <typename FormatContext>
  auto format(std::complex<double> const &number, FormatContext &ctx) {
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
