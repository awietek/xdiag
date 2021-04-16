#include "combinations.h"

namespace hydra {

template <class bit_t>
Combinations<bit_t>::Combinations(int n, int k)
  : n_(n), k_(k),
    size_(combinatorics::binomial(n, k))
{  
  if (k>n)
    HydraLog.err("Error constructing Combinations: k>n");
  else if (k<0)
    HydraLog.err("Error constructing Combinations: k<0");
  else if (n<0)
    HydraLog.err("Error constructing Combinations: n<0");
  else
    {
      bit_t begin = (((bit_t)1 << k) - 1);
      bit_t end = begin << (n - k);
      end = combinatorics::get_next_pattern(end);
      begin_ = CombinationsIterator(begin);
      end_ = CombinationsIterator(end);
    }
}

template <class bit_t>
CombinationsIterator<bit_t>::CombinationsIterator(bit_t state)
    : current_(state) {}

template class Combinations<uint16>;
template class Combinations<uint32>;
template class Combinations<uint64>;

template class CombinationsIterator<uint16>;
template class CombinationsIterator<uint32>;
template class CombinationsIterator<uint64>;

} // namespace hydra
