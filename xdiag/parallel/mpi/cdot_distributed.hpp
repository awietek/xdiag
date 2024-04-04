#pragma once
#ifdef XDIAG_USE_MPI

#include <mpi.h>

#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/common.hpp>

namespace xdiag {

template <class coeff_t>
coeff_t cdot_distributed(arma::Col<coeff_t> const &v,
			 arma::Col<coeff_t> const &w);

template <class coeff_t> double norm_distributed(arma::Col<coeff_t> const &v);

} // namespace xdiag

namespace xdiag::mpi {
using scomplex = std::complex<float>;

double stable_dot_product(uint64_t n, const double *x, const double *y);
double stable_dot_product(uint64_t n, const float *x, const float *y);
complex stable_dot_product(uint64_t n, const complex *x, const complex *y);
scomplex stable_dot_product(uint64_t n, const scomplex *x, const scomplex *y);

} // namespace xdiag::mpi
#endif
