#include <algorithm>
#include <cassert>

#include "index_search.h"

#include <hydra/bases/basis_spinhalf.h>

namespace hydra {

template <class basis_t, class idx_t>
IndexSearch<basis_t, idx_t>::IndexSearch(const basis_t &basis)
    : size_(basis.size()) {
  states_.clear();
  state_t previous;
  int ctr = 0;
  for (auto state : basis) {
    if (ctr++ != 0)
      assert(state > previous);
    previous = state;
    states_.push_back(state);
  }
}

template <class basis_t, class idx_t>
idx_t IndexSearch<basis_t, idx_t>::index(const state_t &state) const {
  // Do binary search
  return std::lower_bound(states_.begin(), states_.end(), state) -
         states_.begin();
}

template class IndexSearch<BasisSpinHalf<uint16>, uint16>;
template class IndexSearch<BasisSpinHalf<uint16>, uint32>;
template class IndexSearch<BasisSpinHalf<uint16>, uint64>;

template class IndexSearch<BasisSpinHalf<uint32>, uint16>;
template class IndexSearch<BasisSpinHalf<uint32>, uint32>;
template class IndexSearch<BasisSpinHalf<uint32>, uint64>;

template class IndexSearch<BasisSpinHalf<uint64>, uint16>;
template class IndexSearch<BasisSpinHalf<uint64>, uint32>;
template class IndexSearch<BasisSpinHalf<uint64>, uint64>;

} // namespace hydra
