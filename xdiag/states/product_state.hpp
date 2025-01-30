#pragma once

#include <string>
#include <vector>
#include <xdiag/common.hpp>

namespace xdiag {

class ProductState {
public:
  using iterator_t = std::vector<std::string>::const_iterator;

  ProductState() = default;
  XDIAG_API explicit ProductState(int64_t nsites);
  XDIAG_API explicit ProductState(std::vector<std::string> const &local_states);

  XDIAG_API std::string const &operator[](int64_t i) const;
  XDIAG_API std::string &operator[](int64_t i);

  int64_t size() const;
  int64_t nsites() const;
  XDIAG_API void push_back(std::string l);

  XDIAG_API iterator_t begin() const;
  XDIAG_API iterator_t end() const;

  XDIAG_API bool operator==(ProductState const &rhs) const;
  XDIAG_API bool operator!=(ProductState const &rhs) const;

private:
  std::vector<std::string> local_states_;
};

XDIAG_API int64_t size(ProductState const &p);
XDIAG_API int64_t nsites(ProductState const &p);
XDIAG_API std::ostream &operator<<(std::ostream &out,
                                   ProductState const &state);
XDIAG_API std::string to_string(ProductState const &state,
                                std::string format = "fancy");

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
