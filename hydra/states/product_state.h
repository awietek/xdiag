#pragma once

#include <string>
#include <vector>

#include <hydra/blocks/blocks.h>
#include <hydra/common.h>
#include <hydra/states/state.h>

namespace hydra {

class ProductState {
public:
  using iterator_t = std::vector<std::string>::const_iterator;

  ProductState() = default;
  explicit ProductState(std::vector<std::string> const &local_states);
  inline std::string const &operator[](int i) const { return local_states_[i]; }
  inline std::string &operator[](int i) { return local_states_[i]; }
  inline int n_sites() const { return local_states_.size(); }
  inline void operator<<(std::string l) { local_states_.push_back(l); }
  inline iterator_t begin() const { return local_states_.begin(); }
  inline iterator_t end() const { return local_states_.end(); }

private:
  std::vector<std::string> local_states_;
};

template <typename coeff_t> class State;

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
