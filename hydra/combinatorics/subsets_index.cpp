#include "subsets_index.h"
#include <lila/utils/logger.h>

namespace hydra::combinatorics {

template <class bit_t>
SubsetsIndex<bit_t>::SubsetsIndex(int n) : n_(n), size_((idx_t)1 << n) {
  if (n < 0)
    lila::Log.err("Error constructing SubsetsIndex: n<0");
  else {
    bit_t begin = 0;
    bit_t end = (idx_t)1 << n;
    begin_ = SubsetsIndexIterator(begin);
    end_ = SubsetsIndexIterator(end);
  }
}

template <class bit_t>
SubsetsIndexIterator<bit_t>::SubsetsIndexIterator(bit_t state) : current_(state) {}

template class SubsetsIndex<uint16_t>;
template class SubsetsIndex<uint32_t>;
template class SubsetsIndex<uint64_t>;

template class SubsetsIndexIterator<uint16_t>;
template class SubsetsIndexIterator<uint32_t>;
template class SubsetsIndexIterator<uint64_t>;

} // namespace hydra
