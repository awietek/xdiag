#include "indexing_no_np.h"

namespace hydra::indexing::electron {

using namespace combinatorics;

template <typename bit_t>
IndexingNoNp<bit_t>::IndexingNoNp(int n_sites)
    : n_sites_(n_sites), size_ups_((idx_t)1 << n_sites),
      size_dns_((idx_t)1 << n_sites), size_(size_ups_ * size_dns_) {}

template <typename bit_t> int IndexingNoNp<bit_t>::n_sites() const {
  return n_sites_;
}
template <typename bit_t> idx_t IndexingNoNp<bit_t>::size_ups() const {
  return size_ups_;
}
template <typename bit_t> idx_t IndexingNoNp<bit_t>::size_dns() const {
  return size_dns_;
}
template <typename bit_t> idx_t IndexingNoNp<bit_t>::size() const {
  return size_;
}

template <typename bit_t>
Subsets<bit_t> IndexingNoNp<bit_t>::states_ups() const {
  return Subsets<bit_t>(n_sites_);
}
template <typename bit_t>
Subsets<bit_t> IndexingNoNp<bit_t>::states_dns() const {
  return Subsets<bit_t>(n_sites_);
}

template <typename bit_t>
SubsetsIndex<bit_t> IndexingNoNp<bit_t>::states_indices_ups() const {
  return SubsetsIndex<bit_t>(n_sites_);
}

template <typename bit_t>
SubsetsIndex<bit_t> IndexingNoNp<bit_t>::states_indices_dns() const {
  return SubsetsIndex<bit_t>(n_sites_);
}

#ifdef _OPENMP
template <typename bit_t>
SubsetsThread<bit_t> IndexingNoNp<bit_t>::states_ups_thread() const {
  return SubsetsThread<bit_t>(n_sites_);
}

template <typename bit_t>
SubsetsThread<bit_t> IndexingNoNp<bit_t>::states_dns_thread() const {
  return SubsetsThread<bit_t>(n_sites_);
}

template <typename bit_t>
SubsetsIndexThread<bit_t> IndexingNoNp<bit_t>::states_indices_ups_thread() const {
  return SubsetsIndexThread<bit_t>(n_sites_);
}

template <typename bit_t>
SubsetsIndexThread<bit_t> IndexingNoNp<bit_t>::states_indices_dns_thread() const {
  return SubsetsIndexThread<bit_t>(n_sites_);
}
#endif

template class IndexingNoNp<uint16_t>;
template class IndexingNoNp<uint32_t>;
template class IndexingNoNp<uint64_t>;

} // namespace hydra::indexing::electron
