#pragma once

#include <string>
#include <vector>
#include <xdiag/common.hpp>

namespace xdiag {

class State;
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

} // namespace xdiag
