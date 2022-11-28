#include "random_state.h"

namespace hydra {

RandomState::RandomState(uint64_t seed) : seed_(seed) {}

template <typename coeff_t>
void fill(RandomState const &rstate, State<coeff_t> &state) {
  uint64_t seed = rstate.seed();

  // the random numbers should be different for different blocks
  uint64_t seed_modified = random::hash_combine(seed, hash(state.block()));
  random::fill_random_normal_vector(state.vector(), seed_modified);

  // normalize
  coeff_t norm = arma::norm(state.vector());
  state.vector() /= norm;
}

template void fill(RandomState const &, State<double> &);
template void fill(RandomState const &, State<complex> &);

template <typename coeff_t>
State<coeff_t> random_state(Block const &block, uint64_t seed) {
  return State<coeff_t>(block, RandomState(seed));
}

template State<double> random_state<double>(Block const &block, uint64_t seed);
template State<complex> random_state<complex>(Block const &block,
                                              uint64_t seed);
StateReal random_state_real(Block const &block, uint64_t seed) {
  return random_state<double>(block, seed);
}

StateCplx random_state_cplx(Block const &block, uint64_t seed) {
  return random_state<complex>(block, seed);
}

} // namespace hydra
