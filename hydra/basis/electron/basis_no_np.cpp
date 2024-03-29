#include "basis_no_np.h"

namespace hydra::basis::electron {

using namespace combinatorics;

template <typename bit_t>
BasisNoNp<bit_t>::BasisNoNp(int n_sites)
    : n_sites_(n_sites), size_ups_((int64_t)1 << n_sites),
      size_dns_((int64_t)1 << n_sites), size_(size_ups_ * size_dns_) {
  if (n_sites < 0) {
    throw(std::invalid_argument("n_sites < 0"));
  }
}

template <typename bit_t> int BasisNoNp<bit_t>::n_sites() const {
  return n_sites_;
}
template <typename bit_t> int64_t BasisNoNp<bit_t>::size_ups() const {
  return size_ups_;
}
template <typename bit_t> int64_t BasisNoNp<bit_t>::size_dns() const {
  return size_dns_;
}
template <typename bit_t> int64_t BasisNoNp<bit_t>::size() const {
  return size_;
}
template <typename bit_t> int64_t BasisNoNp<bit_t>::dim() const {
  return size_;
}

template <typename bit_t> Subsets<bit_t> BasisNoNp<bit_t>::states_ups() const {
  return Subsets<bit_t>(n_sites_);
}
template <typename bit_t> Subsets<bit_t> BasisNoNp<bit_t>::states_dns() const {
  return Subsets<bit_t>(n_sites_);
}

template <typename bit_t>
SubsetsIndex<bit_t> BasisNoNp<bit_t>::states_indices_ups() const {
  return SubsetsIndex<bit_t>(n_sites_);
}

template <typename bit_t>
SubsetsIndex<bit_t> BasisNoNp<bit_t>::states_indices_dns() const {
  return SubsetsIndex<bit_t>(n_sites_);
}

#ifdef _OPENMP
template <typename bit_t>
SubsetsThread<bit_t> BasisNoNp<bit_t>::states_ups_thread() const {
  return SubsetsThread<bit_t>(n_sites_);
}

template <typename bit_t>
SubsetsThread<bit_t> BasisNoNp<bit_t>::states_dns_thread() const {
  return SubsetsThread<bit_t>(n_sites_);
}

template <typename bit_t>
SubsetsIndexThread<bit_t> BasisNoNp<bit_t>::states_indices_ups_thread() const {
  return SubsetsIndexThread<bit_t>(n_sites_);
}

template <typename bit_t>
SubsetsIndexThread<bit_t> BasisNoNp<bit_t>::states_indices_dns_thread() const {
  return SubsetsIndexThread<bit_t>(n_sites_);
}
#endif

template class BasisNoNp<uint16_t>;
template class BasisNoNp<uint32_t>;
template class BasisNoNp<uint64_t>;

} // namespace hydra::basis::electron
