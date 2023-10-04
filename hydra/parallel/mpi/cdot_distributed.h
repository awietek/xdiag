#pragma once
#ifdef HYDRA_USE_MPI

#include "extern/armadillo/armadillo"
#include <hydra/common.h>
#include <mpi.h>

namespace hydra {

template <class coeff_t>
coeff_t cdot_distributed(arma::Col<coeff_t> const &v,
			 arma::Col<coeff_t> const &w);

template <class coeff_t> double norm_distributed(arma::Col<coeff_t> const &v);

} // namespace hydra

namespace hydra::mpi {
using scomplex = std::complex<float>;

double stable_dot_product(uint64_t n, const double *x, const double *y);
double stable_dot_product(uint64_t n, const float *x, const float *y);
complex stable_dot_product(uint64_t n, const complex *x, const complex *y);
scomplex stable_dot_product(uint64_t n, const scomplex *x, const scomplex *y);

} // namespace hydra::mpi
#endif
