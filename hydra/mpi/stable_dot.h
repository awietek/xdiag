#pragma once

#include <lila/all.h>
#include <hydra/common.h>

namespace hydra::mpi {

template <class coeff_t>
coeff_t StableDot(lila::Vector<coeff_t> const &v, lila::Vector<coeff_t> const &w);

template <class coeff_t>
coeff_t StableNorm(lila::Vector<coeff_t> const &v);


namespace detail {

double stable_dot_product(const uint64 &n, const double *x, const double *y);

double stable_dot_product(const uint64 &n, const float *x, const float *y);

complex stable_dot_product(const uint64 &n, const complex *x, const complex *y);

scomplex stable_dot_product(const uint64 &n, const scomplex *x,
                            const scomplex *y);

} // namespace detail
} // namespace hydra::mpi
