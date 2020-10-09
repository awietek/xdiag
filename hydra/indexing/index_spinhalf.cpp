#include "index_spinhalf.h"

#include <hydra/utils/combinatorics.h>

namespace hydra {

template <class bit_t, class index_t>
IndexSpinHalf<bit_t, index_t>::IndexSpinHalf(const basis_t &basis)
    : size_(basis.size()), n_sites_(basis.n_sites()),
      qn_(basis.qn()) {}

template <class bit_t, class index_t>
index_t IndexSpinHalf<bit_t, index_t>::index(const state_t &state) const {
  return combinatorics::get_n_for_pattern((uint64)state.spins, n_sites_,
                                          qn_.n_up);
}

template <class bit_t, class index_t>
state_spinhalf<bit_t>
IndexSpinHalf<bit_t, index_t>::state(const index_t &index) const {
  return state_t({(bit_t)combinatorics::get_nth_pattern((uint64)index, n_sites_,
                                                        qn_.n_up)});
}

template class IndexSpinHalf<uint16, uint16>;
template class IndexSpinHalf<uint16, uint32>;
template class IndexSpinHalf<uint16, uint64>;

template class IndexSpinHalf<uint32, uint16>;
template class IndexSpinHalf<uint32, uint32>;
template class IndexSpinHalf<uint32, uint64>;

template class IndexSpinHalf<uint64, uint16>;
template class IndexSpinHalf<uint64, uint32>;
template class IndexSpinHalf<uint64, uint64>;

} // namespace hydra
