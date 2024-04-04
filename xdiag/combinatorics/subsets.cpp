#include "subsets.hpp"
#include <xdiag/utils/logger.hpp>

#ifdef _OPENMP
#include <xdiag/parallel/omp/omp_utils.hpp>
#endif

namespace xdiag::combinatorics {

template <class bit_t>
Subsets<bit_t>::Subsets(int64_t n) : n_(n), size_((int64_t)1 << n) {
  if (n < 0)
    Log.err("Error constructing Subsets: n<0");
  else {
    begin_ = SubsetsIterator<bit_t>((int64_t)0);
    end_ = SubsetsIterator<bit_t>((int64_t)1 << n);
  }
}

template class Subsets<uint16_t>;
template class Subsets<uint32_t>;
template class Subsets<uint64_t>;

template <class bit_t>
SubsetsIterator<bit_t>::SubsetsIterator(int64_t idx) : current_((bit_t)idx) {}

template class SubsetsIterator<uint16_t>;
template class SubsetsIterator<uint32_t>;
template class SubsetsIterator<uint64_t>;

#ifdef _OPENMP
template <class bit_t>
SubsetsThread<bit_t>::SubsetsThread(int64_t n)
    : n_(n)
// , size_((int64_t)1 << n)
{
  if (n < 0)
    Log.err("Error constructing Subsets: n<0");
  else {
    auto [begin_idx, end_idx] = omp::get_omp_subsets_start_end(n);
    begin_ = SubsetsIterator<bit_t>(begin_idx);
    end_ = SubsetsIterator<bit_t>(end_idx);
    size_ = end_idx - begin_idx;
  }
}
template class SubsetsThread<uint16_t>;
template class SubsetsThread<uint32_t>;
template class SubsetsThread<uint64_t>;
#endif
} // namespace xdiag::combinatorics
