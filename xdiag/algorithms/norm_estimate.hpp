#pragma once

#include <functional>

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag {

// Returns an estimate of the 1-norm of an operator (needed by some iterative
// algorithms)

double norm_estimate(OpSum const &ops, Block const &block,
                     int64_t n_max_attempts = 5, uint64_t seed = 42);

template <typename block_t>
double norm_estimate(OpSum const &ops, block_t const &block,
                     int64_t n_max_attempts = 5, uint64_t seed = 42);

template <typename coeff_t>
double norm_estimate(arma::Mat<coeff_t> const &A, int64_t n_max_attempts = 5,
                     uint64_t seed = 42); // used mainly for testing

} // namespace xdiag
