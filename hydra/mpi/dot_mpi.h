#pragma once

#include <hydra/common.h>
#include <lila/all.h>

namespace hydra {

template <class coeff_t>
coeff_t DotMPI(lila::Vector<coeff_t> const &v, lila::Vector<coeff_t> const &w);

template <class coeff_t>
lila::real_t<coeff_t> NormMPI(lila::Vector<coeff_t> const &v);

} // namespace hydra

namespace hydra::mpi {

double stable_dot_product(uint64_t n, const double *x, const double *y);
double stable_dot_product(uint64_t n, const float *x, const float *y);
complex stable_dot_product(uint64_t n, const complex *x, const complex *y);
scomplex stable_dot_product(uint64_t n, const scomplex *x, const scomplex *y);

} // namespace hydra::mpi
