// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/armadillo.hpp>
#include <xdiag/math/complex.hpp>

namespace xdiag {

bool isapprox(double a, double b, double rtol = 1e-12, double atol = 1e-12);
bool isapprox(complex a, complex b, double rtol = 1e-12, double atol = 1e-12);

template <typename coeff_t>
bool isapprox(arma::Col<coeff_t> const &v, arma::Col<coeff_t> const &w,
              double rtol = 1e-12, double atol = 1e-12);

template <typename coeff_t>
bool isapprox(arma::Mat<coeff_t> const &v, arma::Mat<coeff_t> const &w,
              double rtol = 1e-12, double atol = 1e-12);

} // namespace xdiag
