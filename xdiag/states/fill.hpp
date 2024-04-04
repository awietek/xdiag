#pragma once

#include <xdiag/common.hpp>
#include <xdiag/states/state.hpp>

namespace xdiag {

template <class fill_f>
void fill(State &state, fill_f fun, int64_t col = 0) try {
  auto const &block = state.block();
  if (col > state.n_cols()) {
    throw(std::invalid_argument(
        "Column number larger that number of colums in state"));
  }

  if (state.isreal()) {
    arma::vec vec = state.vector(col, false);
    std::visit([fun, &vec](auto &&block) { fill(block, fun, vec); }, block);
  } else {
    arma::cx_vec vec = state.vectorC(col, false);
    std::visit([fun, &vec](auto &&block) { fill(block, fun, vec); }, block);
  }
} catch (...) {
  rethrow(__func__, std::runtime_error("Unable to fill State"));
}

} // namespace xdiag
