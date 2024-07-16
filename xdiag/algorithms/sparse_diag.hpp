#pragma once

#include <tuple>

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/state.hpp>

namespace xdiag {

double eigval0(OpSum const &ops, block_variant_t const &block,
               double precision = 1e-12, int64_t max_iterations = 1000,
               bool force_complex = false, int64_t random_seed = 42);

std::tuple<double, State>
eig0(OpSum const &ops, block_variant_t const &block,
     double precision = 1e-12, int64_t max_iterations = 1000,
     bool force_complex = false, int64_t random_seed = 42);

} // namespace xdiag
