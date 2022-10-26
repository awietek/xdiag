#pragma once

#include <string>
#include <vector>

#include <hydra/blocks/blocks.h>
#include <hydra/common.h>
#include <hydra/states/state.h>

namespace hydra {

class ProductState {
public:
  ProductState() = default;
  explicit ProductState(std::vector<std::string> const &local_states);
  inline std::string operator[](int i) const { return local_states_[i]; }
  inline int n_sites() const { return local_states_.size(); }

private:
  std::vector<std::string> local_states_;
};

template <typename coeff_t>
void fill(ProductState const &pstate, State<coeff_t> &state);

template <typename coeff_t>
void fill(ProductState const &pstate, Spinhalf const &block,
          arma::Col<coeff_t> &vector);
template <typename coeff_t>
void fill(ProductState const &pstate, tJ const &block,
          arma::Col<coeff_t> &vector);
template <typename coeff_t>
void fill(ProductState const &pstate, Electron const &block,
          arma::Col<coeff_t> &vector);

} // namespace hydra
