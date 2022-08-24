#include "subsets.h"
#include <hydra/utils/logger.h>

#ifdef HYDRA_ENABLE_OPENMP
#include <hydra/utils/openmp_utils.h>
#endif

namespace hydra::combinatorics {

template <class bit_t>
Subsets<bit_t>::Subsets(int n) : n_(n), size_((idx_t)1 << n) {
  if (n < 0)
    Log.err("Error constructing Subsets: n<0");
  else {
    begin_ = SubsetsIterator<bit_t>((idx_t)0);
    end_ = SubsetsIterator<bit_t>((idx_t)1 << n);
  }
}

template class Subsets<uint16_t>;
template class Subsets<uint32_t>;
template class Subsets<uint64_t>;

template <class bit_t>
SubsetsIterator<bit_t>::SubsetsIterator(idx_t idx) : current_((bit_t)idx) {}

template class SubsetsIterator<uint16_t>;
template class SubsetsIterator<uint32_t>;
template class SubsetsIterator<uint64_t>;

#ifdef HYDRA_ENABLE_OPENMP
template <class bit_t>
SubsetsThread<bit_t>::SubsetsThread(int n)
    : n_(n)
// , size_((idx_t)1 << n)
{
  if (n < 0)
    Log.err("Error constructing Subsets: n<0");
  else {
    auto [begin_idx, end_idx] = utils::get_omp_subsets_start_end(n);
    begin_ = SubsetsIterator<bit_t>(begin_idx);
    end_ = SubsetsIterator<bit_t>(end_idx);
    size_ = end_idx - begin_idx;
  }
}
template class SubsetsThread<uint16_t>;
template class SubsetsThread<uint32_t>;
template class SubsetsThread<uint64_t>;
#endif
} // namespace hydra::combinatorics
