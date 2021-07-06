#pragma once

#include <complex>
#include <lila/all.h>

namespace hydra {

using int16 = short;
using int32 = int;
using int64 = long;
using uint16 = unsigned short;
using uint32 = unsigned int;
using uint64 = unsigned long;

using std_bit_t = uint64;
using number_t = int32;

using idx_t = int64;

using scomplex = std::complex<float>;
using complex = std::complex<double>;

template <class T> struct is_complex_t : public std::false_type {};
template <class T>
struct is_complex_t<std::complex<T>> : public std::true_type {};
template <class T> constexpr bool is_complex() {
  return is_complex_t<T>::value;
}

constexpr idx_t invalid_index = -1;
constexpr bool index_not_found(idx_t idx) { return idx < 0; }
constexpr bool index_valid(idx_t idx) { return idx >= 0; }

} // namespace hydra
