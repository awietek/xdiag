#pragma once

#include <xdiag/common.hpp>
#include <xdiag/states/gpwf.hpp>
#include <xdiag/states/product_state.hpp>
#include <xdiag/states/random_state.hpp>
#include <xdiag/states/state.hpp>

namespace xdiag {

void fill(State &state, std::function<double(ProductState const &)> coeff_f,
          int64_t col = 0);
void fill(State &state, std::function<complex(ProductState const &)> coeff_f,
          int64_t col = 0);

void fill(State &state, RandomState const &rstate, int64_t col = 0);
void fill(State &state, ProductState const &pstate, int64_t col = 0);
void fill(State &state, GPWF const &gpwf, int64_t col = 0);

} // namespace xdiag
