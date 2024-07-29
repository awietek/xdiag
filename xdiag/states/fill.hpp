#pragma once

#include <xdiag/common.hpp>
#include <xdiag/states/product_state.hpp>
#include <xdiag/states/random_state.hpp>
#include <xdiag/states/state.hpp>

namespace xdiag {

void fill(State &state, RandomState const &rstate, int64_t col = 0);
void fill(State &state, ProductState const &pstate, int64_t col = 0);
template <class block_t, typename coeff_t>
void fill(block_t const &block, arma::Col<coeff_t> &vec,
          ProductState const &pstate);

} // namespace xdiag
