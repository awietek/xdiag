#include "index_table.h"

#include <hydra/bases/basis_spinhalf.h>

namespace hydra {

template <class hilbertspace_t, class idx_t>
IndexTable<hilbertspace_t, idx_t>::IndexTable(
    const hilbertspace_t &hilbertspace)
    : size_(hilbertspace.size()), n_sites_(hilbertspace.n_sites()),
      indices_(hilbertspace.rawsize(), -1) {
  states_.clear();
  idx_t idx = 0;

  for (auto state : hilbertspace) {
    states_.push_back(state);
    indices_[rawidx(state, n_sites_)] = idx;
    ++idx;
  }
}

template class IndexTable<BasisSpinHalf<uint16>, uint16>;
template class IndexTable<BasisSpinHalf<uint16>, uint32>;
template class IndexTable<BasisSpinHalf<uint16>, uint64>;

template class IndexTable<BasisSpinHalf<uint32>, uint16>;
template class IndexTable<BasisSpinHalf<uint32>, uint32>;
template class IndexTable<BasisSpinHalf<uint32>, uint64>;

template class IndexTable<BasisSpinHalf<uint64>, uint16>;
template class IndexTable<BasisSpinHalf<uint64>, uint32>;
template class IndexTable<BasisSpinHalf<uint64>, uint64>;

} // namespace hydra
