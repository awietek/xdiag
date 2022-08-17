#include "subsets_index.h"

#include <hydra/utils/logger.h>

#ifdef HYDRA_ENABLE_OPENMP
#include <hydra/utils/openmp_utils.h>
#endif

namespace hydra::combinatorics {

template <class bit_t>
SubsetsIndex<bit_t>::SubsetsIndex(int n) : n_(n), size_((idx_t)1 << n) {
  if (n < 0) {
    Log.err("Error constructing SubsetsIndex: n<0");
  } else {
    begin_ = SubsetsIndexIterator<bit_t>(0);
    end_ = SubsetsIndexIterator<bit_t>((idx_t)1 << n);
  }
}

template class SubsetsIndex<uint16_t>;
template class SubsetsIndex<uint32_t>;
template class SubsetsIndex<uint64_t>;

template <class bit_t>
SubsetsIndexIterator<bit_t>::SubsetsIndexIterator(idx_t idx)
    : current_((bit_t)idx) {}

template class SubsetsIndexIterator<uint16_t>;
template class SubsetsIndexIterator<uint32_t>;
template class SubsetsIndexIterator<uint64_t>;

#ifdef HYDRA_ENABLE_OPENMP
template <class bit_t>
SubsetsIndexThread<bit_t>::SubsetsIndexThread(int n)
    : n_(n), size_((idx_t)1 << n) {
  if (n < 0) {
    Log.err("Error constructing SubsetsIndex: n<0");
  } else {
    auto [begin_idx, end_idx] = get_omp_subsets_start_end(n);
    begin_ = SubsetsIndexIterator<bit_t>(begin_idx);
    end_ = SubsetsIndexIterator<bit_t>(end_idx);
  }
}

template class SubsetsIndexThread<uint16_t>;
template class SubsetsIndexThread<uint32_t>;
template class SubsetsIndexThread<uint64_t>;
#endif

} // namespace hydra::combinatorics
