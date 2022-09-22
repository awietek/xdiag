#pragma once

#include <complex>
#include <cstdint>
#include <utility>

#include <hydra/bitops/bitops.h>
#define BSTR(x) bitops::bits_to_string(x, n_sites)

#include <hydra/utils/logger.h>

#include "extern/armadillo/armadillo"

namespace hydra {

using std_bit_t = uint64_t;
using number_t = int32_t;

using idx_t = int64_t;

using scomplex = std::complex<float>;
using complex = std::complex<double>;

template <class T> struct is_complex_t : public std::false_type {};
template <class T>
struct is_complex_t<std::complex<T>> : public std::true_type {};
template <class T> constexpr bool is_complex() {
  return is_complex_t<T>::value;
}
template <class T> constexpr bool is_real() { return !is_complex_t<T>::value; }

constexpr idx_t invalid_index = (idx_t)-1;
constexpr int invalid_n = (idx_t)-1;

constexpr int undefined_qn = std::numeric_limits<int>::min();
constexpr std::pair<int, int> undefined_qns = {undefined_qn, undefined_qn};

constexpr bool index_not_found(idx_t idx) { return idx < 0; }
constexpr bool index_valid(idx_t idx) { return idx >= 0; }

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

template <class coeff_t>
using real_t = typename real_type_struct<coeff_t>::type;

template <class coeff_t>
using complex_t = typename complex_type_struct<coeff_t>::type;

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

} // namespace hydra
