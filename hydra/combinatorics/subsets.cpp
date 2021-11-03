#include "subsets.h"
#include <lila/utils/logger.h>

namespace hydra::combinatorics {

template <class bit_t>
Subsets<bit_t>::Subsets(int n) : n_(n), size_((idx_t)1 << n) {
  if (n < 0)
    lila::Log.err("Error constructing Subsets: n<0");
  else {
    bit_t begin = 0;
    bit_t end = (idx_t)1 << n;
    begin_ = SubsetsIterator(begin);
    end_ = SubsetsIterator(end);
  }
}

template <class bit_t>
SubsetsIterator<bit_t>::SubsetsIterator(bit_t state) : current_(state) {}

template class Subsets<uint16_t>;
template class Subsets<uint32_t>;
template class Subsets<uint64_t>;

template class SubsetsIterator<uint16_t>;
template class SubsetsIterator<uint32_t>;
template class SubsetsIterator<uint64_t>;

} // namespace hydra
