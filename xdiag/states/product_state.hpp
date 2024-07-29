#pragma once

#include <string>
#include <vector>
#include <xdiag/common.hpp>

namespace xdiag {

class ProductState {
public:
  using iterator_t = std::vector<std::string>::const_iterator;

  ProductState() = default;
  explicit ProductState(int64_t n_sites);
  explicit ProductState(std::vector<std::string> const &local_states);

  std::string const &operator[](int64_t i) const;
  std::string &operator[](int64_t i);

  int64_t n_sites() const;
  void operator<<(std::string l);
  iterator_t begin() const;
  iterator_t end() const;

  bool operator==(ProductState const &rhs) const;
  bool operator!=(ProductState const &rhs) const;

private:
  std::vector<std::string> local_states_;
};

std::ostream &operator<<(std::ostream &out, ProductState const &state);
std::string to_string(ProductState const &state, std::string format = "fancy");

// Conversions for low-level encoding of states
template <typename bit_t>
void to_product_state_spinhalf(bit_t spins, ProductState &pstate);
template <typename bit_t>
void to_product_state_tj(bit_t ups, bit_t dns, ProductState &pstate);
template <typename bit_t>
void to_product_state_electron(bit_t ups, bit_t dns, ProductState &pstate);

template <typename bit_t> bit_t to_bits_spinhalf(ProductState const &pstate);
template <typename bit_t>
std::pair<bit_t, bit_t> to_bits_tj(ProductState const &pstate);
template <typename bit_t>
std::pair<bit_t, bit_t> to_bits_electron(ProductState const &pstate);

} // namespace xdiag
