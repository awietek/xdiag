#include "create_state.hpp"

#include <xdiag/states/fill.hpp>

namespace xdiag {
State product(Block const &block, std::vector<std::string> const &local_state,
              bool real) {
  return std::visit(
      [&](auto &&block) { return product(block, local_state, real); }, block);
}

template <typename block_t>
State product(block_t const &block, std::vector<std::string> const &local_state,
              bool real) {
  auto state = State(block, real);
  auto pstate = ProductState(local_state);
  fill(state, pstate);
  return state;
}

template State product(Spinhalf const &, std::vector<std::string> const &,
                       bool);
template State product(tJ const &, std::vector<std::string> const &, bool);
template State product(Electron const &, std::vector<std::string> const &,
                       bool);

#ifdef XDIAG_USE_MPI
template State product(SpinhalfDistributed const &,
                       std::vector<std::string> const &, bool);
template State product(tJDistributed const &, std::vector<std::string> const &,
                       bool);
#endif

State rand(Block const &block, bool real, int64_t seed, bool normalized) {
  return std::visit(
      [&](auto &&block) { return rand(block, real, seed, normalized); }, block);
}

template <typename block_t>
State rand(block_t const &block, bool real, int64_t seed, bool normalized) {
  auto state = State(block, real);
  auto rstate = RandomState(seed, normalized);
  fill(state, rstate);
  return state;
}
template State rand(Spinhalf const &, bool, int64_t, bool);
template State rand(tJ const &, bool, int64_t, bool);
template State rand(Electron const &, bool, int64_t, bool);
#ifdef XDIAG_USE_MPI
template State rand(SpinhalfDistributed const &, bool, int64_t, bool);
template State rand(tJDistributed const &, bool, int64_t, bool);
#endif

State zeros(Block const &block, bool real, int64_t n_cols) {
  return std::visit([&](auto &&blk) { return zeros(blk, real, n_cols); },
                    block);
}

template <typename block_t>
State zeros(block_t const &block, bool real, int64_t n_cols) {
  return State(block, real, n_cols);
}
template State zeros(Spinhalf const &, bool, int64_t);
template State zeros(tJ const &, bool, int64_t);
template State zeros(Electron const &, bool, int64_t);
#ifdef XDIAG_USE_MPI
template State zeros(tJDistributed const &, bool, int64_t);
template State zeros(SpinhalfDistributed const &, bool, int64_t);
#endif

void zero(State &state) {
  if (state.isreal()) {
    state.matrix(false).zeros();
  } else {
    state.matrixC(false).zeros();
  }
}

} // namespace xdiag
