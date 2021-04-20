#include "subsets.h"

namespace hydra {
namespace combinatorics {

template <class bit_t>
Subsets<bit_t>::Subsets(int n)
  : n_(n),
    size_((idx_t)1 << n)
{
  if (n<0)
    HydraLog.err("Error constructing Subsets: n<0");
  else
    {
      bit_t begin = 0;
      bit_t end = (idx_t)1 << n;
      begin_ = SubsetsIterator(begin);
      end_ = SubsetsIterator(end);
    }
}

template <class bit_t>
SubsetsIterator<bit_t>::SubsetsIterator(bit_t state)
    : current_(state) {}

template class Subsets<uint16>;
template class Subsets<uint32>;
template class Subsets<uint64>;

template class SubsetsIterator<uint16>;
template class SubsetsIterator<uint32>;
template class SubsetsIterator<uint64>;

} // namespace combinatorics
} // namespace hydra
