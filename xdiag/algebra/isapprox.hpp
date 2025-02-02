#pragma once

#include <xdiag/common.hpp>
#include <xdiag/extern/armadillo/armadillo>

namespace xdiag {

XDIAG_API bool isapprox(double x, double y, double rtol = 1e-12,
                        double atol = 1e-12);
XDIAG_API bool isapprox(complex x, complex y, double rtol = 1e-12,
                        double atol = 1e-12);

XDIAG_API bool isapprox(arma::vec const &v, arma::vec const &w,
                        double rtol = 1e-12, double atol = 1e-12);
XDIAG_API bool isapprox(arma::cx_vec const &v, arma::cx_vec const &w,
                        double rtol = 1e-12, double atol = 1e-12);

XDIAG_API bool isapprox(arma::mat const &A, arma::mat const &B,
                        double rtol = 1e-12, double atol = 1e-12);
XDIAG_API bool isapprox(arma::cx_mat const &A, arma::cx_mat const &B,
                        double rtol = 1e-12, double atol = 1e-12);

} // namespace xdiag
