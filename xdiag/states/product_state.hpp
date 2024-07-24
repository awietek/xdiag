#pragma once

#include <string>
#include <vector>

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/common.hpp>
#include <xdiag/states/state.hpp>

namespace xdiag {

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

State product_state(Block const &block,
                    std::vector<std::string> const &local_state,
                    bool real = true);

template <typename block_t>
State product_state(block_t const &block,
                    std::vector<std::string> const &local_state,
                    bool real = true);

void fill(State &state, ProductState const &pstate, int64_t col = 0);

template <typename coeff_t>
void fill(Spinhalf const &block, arma::Col<coeff_t> &vec,
          ProductState const &pstate);

template <typename coeff_t>
void fill(tJ const &block, arma::Col<coeff_t> &vec, ProductState const &pstate);

template <typename coeff_t>
void fill(Electron const &block, arma::Col<coeff_t> &vec,
          ProductState const &pstate);

#ifdef XDIAG_USE_MPI
template <typename coeff_t>
void fill(tJDistributed const &block, arma::Col<coeff_t> &vec,
          ProductState const &pstate);
#endif

std::ostream &operator<<(std::ostream &out, ProductState const &state);
std::string to_string(ProductState const &state);
  
} // namespace xdiag
