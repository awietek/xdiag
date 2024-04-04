#pragma once

#include <functional>

#include <xdiag/blocks/blocks.h>
#include <xdiag/common.h>
#include <xdiag/operators/bondlist.h>

namespace xdiag {

// Returns an estimate of the 1-norm of an operator (needed by some iterative
// algorithms)

double norm_estimate(BondList const &bonds, block_variant_t const &block,
                     int64_t n_max_attempts = 5, uint64_t seed = 42);

template <typename block_t>
double norm_estimate(BondList const &bonds, block_t const &block,
                     int64_t n_max_attempts = 5, uint64_t seed = 42);

template <typename coeff_t>
double norm_estimate(arma::Mat<coeff_t> const &A, int64_t n_max_attempts = 5,
                     uint64_t seed = 42); // used mainly for testing

} // namespace xdiag
