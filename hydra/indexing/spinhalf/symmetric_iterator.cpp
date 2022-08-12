#include "symmetric_iterator.h"

namespace hydra::indexing::spinhalf {

template <typename bit_t>
SymmetricIterator<bit_t>::SymmetricIterator(std::vector<bit_t> const &reps,
                                            idx_t idx)
    : data_(reps.data() + idx), idx_(idx) {}

template class SymmetricIterator<uint16_t>;
template class SymmetricIterator<uint32_t>;
template class SymmetricIterator<uint64_t>;

} // namespace hydra::indexing::spinhalf
