// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <complex>

namespace xdiag {

// The two scalar types used throughout xdiag: double (real) and complex.
// All numeric types in Scalar, Vector, and Matrix are variants over these two.
using complex = std::complex<double>;

// Compile-time predicate: isreal<double>() == true, isreal<complex>() == false.
template <typename T> constexpr bool isreal() {
  return std::is_same<double, T>::value;
}

// Uniform real/imag/conj overloads for both double and complex.
// For double: imag returns 0, conj is identity.
constexpr double real(double x) { return x; }
constexpr double real(complex x) { return x.real(); }

constexpr double imag(double) { return 0.; }
constexpr double imag(complex x) { return x.imag(); }

constexpr double conj(double x) { return x; }
constexpr complex conj(complex x) { return complex(x.real(), -x.imag()); }

} // namespace xdiag
