#include "random_state.hpp"

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/random/hash.hpp>
#include <xdiag/random/hash_functions.hpp>
#include <xdiag/random/random_utils.hpp>

namespace xdiag {

RandomState::RandomState(int64_t seed, bool normalized)
    : seed_(seed), normalized_(normalized) {}
int64_t RandomState::seed() const { return seed_; }
bool RandomState::normalized() const { return normalized_; }

void fill(State &state, RandomState const &rstate, int64_t col) try {
  int64_t seed = rstate.seed();
  int64_t seed_modified =
      random::hash_combine(seed, random::hash(state.block()));

  if (state.isreal()) {
    auto v = state.vector(col, false);
    random::fill_random_normal_vector(v, seed_modified);
  } else {
    auto v = state.vectorC(col, false);
    random::fill_random_normal_vector(v, seed_modified);
  }
  if (rstate.normalized()) {
    double nrm = norm(state);
    state /= nrm;
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

State random_state(block_variant_t const &block, bool real, int64_t seed,
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
template State random_state(tJDistributed const &, bool, int64_t, bool);
#endif

} // namespace xdiag
