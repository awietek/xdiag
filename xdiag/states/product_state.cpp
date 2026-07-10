// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "product_state.hpp"

#include <ostream>
#include <sstream>
#include <xdiag/utils/to_string_generic.hpp>

namespace xdiag {

ProductState::ProductState(int64_t nsites) : local_states_(nsites) {}

ProductState::ProductState(std::vector<int32_t> const &local_states)
    : ProductState(
          std::vector<int64_t>(local_states.begin(), local_states.end())) {}

ProductState::ProductState(std::vector<int64_t> const &local_states)
    : local_states_(local_states) {}
ProductState::ProductState(arma::Col<int64_t> const &local_states)
    : ProductState(
          std::vector<int64_t>(local_states.begin(), local_states.end())) {}

int64_t ProductState::operator[](int64_t i) const { return local_states_[i]; }
int64_t &ProductState::operator[](int64_t i) { return local_states_[i]; }

void ProductState::push_back(int64_t l) { local_states_.push_back(l); }
int64_t ProductState::size() const { return local_states_.size(); }
int64_t ProductState::nsites() const { return local_states_.size(); }

ProductState::iterator_t ProductState::begin() const {
  return local_states_.begin();
}
ProductState::iterator_t ProductState::end() const {
  return local_states_.end();
}

bool ProductState::operator==(ProductState const &rhs) const {
  return local_states_ == rhs.local_states_;
}
bool ProductState::operator!=(ProductState const &rhs) const {
  return !operator==(rhs);
}

int64_t size(ProductState const &p) { return p.size(); }
int64_t nsites(ProductState const &p) { return p.nsites(); }

std::ostream &operator<<(std::ostream &out, ProductState const &state) {
  for (int64_t i = 0; i < state.size(); ++i) {
    out << state[i];
    if (i < state.size() - 1) {
      out << " ";
    }
  }
  return out;
}
std::string to_string(ProductState const &state) {
  return utils::to_string_generic(state);
}

} // namespace xdiag
