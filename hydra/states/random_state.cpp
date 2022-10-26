#include "random_state.h"

namespace hydra {

RandomState::RandomState(uint32_t seed) : seed_(seed) {}

template <typename coeff_t>
void fill(RandomState const &rstate, State<coeff_t> &state) {
  uint32_t seed = rstate.seed();

  // the random numbers should be different for different blocks
  uint32_t seed_modified = random::hash_combine(seed, hash(state.block()));
  random::fill_random_normal_vector(state.vector(), seed_modified);

  // normalize
  coeff_t norm = arma::norm(state.vector());
  state.vector() /= norm;
}

} // namespace hydra
