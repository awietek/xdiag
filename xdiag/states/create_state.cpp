#include "create_state.hpp"

#include <xdiag/states/fill.hpp>

namespace xdiag {
State product_state(Block const &block,
                    std::vector<std::string> const &local_state, bool real) {
  return std::visit(
      [&](auto &&block) { return product_state(block, local_state, real); },
      block);
}

template <typename block_t>
State product_state(block_t const &block,
                    std::vector<std::string> const &local_state, bool real) {
  auto state = State(block, real);
  auto pstate = ProductState(local_state);
  fill(state, pstate);
  return state;
}

template State product_state(Spinhalf const &, std::vector<std::string> const &,
                             bool);
template State product_state(tJ const &, std::vector<std::string> const &,
                             bool);
template State product_state(Electron const &, std::vector<std::string> const &,
                             bool);

#ifdef XDIAG_USE_MPI
template State product_state(SpinhalfDistributed const &,
                             std::vector<std::string> const &, bool);
template State product_state(tJDistributed const &,
                             std::vector<std::string> const &, bool);
#endif

State random_state(Block const &block, bool real, int64_t seed,
                   bool normalized) {
  return std::visit(
      [&](auto &&block) { return random_state(block, real, seed, normalized); },
      block);
}

template <typename block_t>
State random_state(block_t const &block, bool real, int64_t seed,
                   bool normalized) {
  auto state = State(block, real);
  auto rstate = RandomState(seed, normalized);
  fill(state, rstate);
  return state;
}
template State random_state(Spinhalf const &, bool, int64_t, bool);
template State random_state(tJ const &, bool, int64_t, bool);
template State random_state(Electron const &, bool, int64_t, bool);
#ifdef XDIAG_USE_MPI
template State random_state(SpinhalfDistributed const &, bool, int64_t, bool);
template State random_state(tJDistributed const &, bool, int64_t, bool);
#endif

State zero_state(Block const &block, bool real, int64_t n_cols) {
  return std::visit([&](auto &&blk) { return zero_state(blk, real, n_cols); },
                    block);
}

template <typename block_t>
State zero_state(block_t const &block, bool real, int64_t n_cols) {
  return State(block, real, n_cols);
}
template State zero_state(Spinhalf const &, bool, int64_t);
template State zero_state(tJ const &, bool, int64_t);
template State zero_state(Electron const &, bool, int64_t);
#ifdef XDIAG_USE_MPI
template State zero_state(tJDistributed const &, bool, int64_t);
template State zero_state(SpinhalfDistributed const &, bool, int64_t);
#endif

void zero(State &state) {
  if (state.isreal()) {
    state.matrix(false).zeros();
  } else {
    state.matrixC(false).zeros();
  }
}

} // namespace xdiag
