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

  std::string const &operator[](int64_t i) const;
  std::string &operator[](int64_t i);

  int64_t n_sites() const;
  void operator<<(std::string l);
  iterator_t begin() const;
  iterator_t end() const;

private:
  std::vector<std::string> local_states_;
};

void fill(State &state, ProductState const &pstate, int64_t col = 0);

template <typename coeff_t>
void fill(Spinhalf const &block, arma::Col<coeff_t> &vec,
          ProductState const &pstate);

template <typename coeff_t>
void fill(tJ const &block, arma::Col<coeff_t> &vec, ProductState const &pstate);

template <typename coeff_t>
void fill(Electron const &block, arma::Col<coeff_t> &vec,
          ProductState const &pstate);

State product_state(block_variant_t const &block,
                    std::vector<std::string> const &local_state,
                    bool real = true);

template <typename block_t>
State product_state(block_t const &block,
                    std::vector<std::string> const &local_state,
                    bool real = true);

} // namespace hydra
