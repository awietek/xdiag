#pragma once

#include <string>
#include <vector>

#include <hydra/blocks/blocks.h>
#include <hydra/common.h>
#include <hydra/states/state.h>
#include <hydra/states/zero_state.h>

namespace hydra {

class ProductState {
public:
  ProductState() = default;
  explicit ProductState(std::vector<std::string> const &local_states);
  inline std::string operator[](int i) { return local_states_[i]; }

private:
  std::vector<std::string> local_states_;
};


void fill(State<coeff_t, Spinhalf> & state, ProductState const& product_state)


  

  
template <typename coeff_t = complex, typename bit_t = std_bit_t>
State<coeff_t, Spinhalf<bit_t>>
product_state(std::vector<std::string> const &local_states,
              Spinhalf<bit_t> const &block);

template <typename bit_t = std_bit_t>
inline StateReal<Spinhalf<bit_t>>
product_state_real(std::vector<std::string> const &local_states,
                   Spinhalf<bit_t> const &block) {
  return product_state<bit_t, double>(local_states, block);
}

template <typename bit_t = std_bit_t>
inline StateCplx<Spinhalf<bit_t>>
product_state_cplx(std::vector<std::string> const &local_states,
                   Spinhalf<bit_t> const &block) {
  return product_state<bit_t, complex>(local_states, block);
}

template <typename coeff_t = complex, typename bit_t = std_bit_t>
State<coeff_t, Electron<bit_t>>
product_state(std::vector<std::string> const &local_states,
              Electron<bit_t> const &block);

template <typename bit_t = std_bit_t>
inline StateReal<Electron<bit_t>>
product_state_real(std::vector<std::string> const &local_states,
                   Electron<bit_t> const &block) {
  return product_state<bit_t, double>(local_states, block);
}

template <typename bit_t = std_bit_t>
inline StateCplx<Electron<bit_t>>
product_state_cplx(std::vector<std::string> const &local_states,
                   Electron<bit_t> const &block) {
  return product_state<bit_t, complex>(local_states, block);
}

template <typename coeff_t = complex, typename bit_t = std_bit_t>
State<coeff_t, tJ<bit_t>>
product_state(std::vector<std::string> const &local_states,
              tJ<bit_t> const &block);

template <typename bit_t = std_bit_t>
inline StateReal<tJ<bit_t>>
product_state_real(std::vector<std::string> const &local_states,
                   tJ<bit_t> const &block) {
  return product_state<bit_t, double>(local_states, block);
}

template <typename bit_t = std_bit_t>
inline StateCplx<tJ<bit_t>>
product_state_cplx(std::vector<std::string> const &local_states,
                   tJ<bit_t> const &block) {
  return product_state<bit_t, complex>(local_states, block);
}
} // namespace hydra
