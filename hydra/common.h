#pragma once

#include <complex>
#include <cstdint>
#include <lila/all.h>
#include <utility>

#include <hydra/bitops/bitops.h>
#define BSTR(x) bitops::bits_to_string(x, n_sites)

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

template <class coeff_t> inline coeff_t complex_to(complex const &cplx) {
  if constexpr (is_complex<coeff_t>())
    return cplx;
  else
    return lila::real(cplx);
}

constexpr idx_t invalid_index = (idx_t)-1;
constexpr int invalid_n = (idx_t)-1;

constexpr int undefined_qn = std::numeric_limits<int>::min();
constexpr std::pair<int, int> undefined_qns = {undefined_qn, undefined_qn};

constexpr bool index_not_found(idx_t idx) { return idx < 0; }
constexpr bool index_valid(idx_t idx) { return idx >= 0; }

// helper type for visitors
template <class... Ts> struct overloaded : Ts... { using Ts::operator()...; };
template <class... Ts> overloaded(Ts...) -> overloaded<Ts...>;
   
} // namespace hydra
