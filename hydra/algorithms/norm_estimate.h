#pragma once

#include <functional>

#include <hydra/blocks/blocks.h>
#include <hydra/common.h>
#include <hydra/operators/bondlist.h>

namespace hydra {

double norm_estimate(BondList const &bonds, block_variant_t const &block);

template <typename block_t>
double norm_estimate(BondList const &bonds, block_t const &block);

double
norm_estimate_real(std::function<arma::vec(arma::vec const &)> const &apply_A,
                   std::function<arma::vec(arma::vec const &)> const &apply_A_T,
                   int64_t N, int n_max_attempts = 5, uint64_t seed = 42);

double norm_estimate_cplx(
    std::function<arma::cx_vec(arma::cx_vec const &)> const &apply_A,
    std::function<arma::cx_vec(arma::cx_vec const &)> const &apply_A_T, int64_t N,
    int n_max_attempts = 5, uint64_t seed = 42);

template <typename coeff_t>
double norm_estimate(
    std::function<arma::Col<coeff_t>(arma::Col<coeff_t> const &)> const
        &apply_A,
    std::function<arma::Col<coeff_t>(arma::Col<coeff_t> const &)> const
        &apply_A_T,
    int64_t N, int n_max_attempts = 5, uint64_t seed = 42);

} // namespace hydra
