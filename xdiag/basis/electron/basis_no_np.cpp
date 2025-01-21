#include "basis_no_np.hpp"

namespace xdiag::basis::electron {

using namespace combinatorics;

template <typename bit_t>
BasisNoNp<bit_t>::BasisNoNp(int n_sites) try
    : n_sites_(n_sites), size_ups_((int64_t)1 << n_sites),
      size_dns_((int64_t)1 << n_sites), size_(size_ups_ * size_dns_) {
  check_n_sites_work_with_bits<bit_t>(n_sites_);
  if (n_sites < 0) {
    XDIAG_THROW("n_sites < 0");
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
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
template <typename bit_t>
typename BasisNoNp<bit_t>::iterator_t BasisNoNp<bit_t>::begin() const {
  return iterator_t(n_sites_, true);
}
template <typename bit_t>
typename BasisNoNp<bit_t>::iterator_t BasisNoNp<bit_t>::end() const {
  return iterator_t(n_sites_, false);
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

template class BasisNoNp<uint32_t>;
template class BasisNoNp<uint64_t>;

template <typename bit_t>
BasisNoNpIterator<bit_t>::BasisNoNpIterator(int64_t n_sites, bool begin)
    : max_((bit_t)1 << n_sites), ups_(begin ? 0 : max_), dns_(0) {}

template <typename bit_t>
BasisNoNpIterator<bit_t> &BasisNoNpIterator<bit_t>::operator++() {
  ++dns_;
  if (dns_ == max_) {
    dns_ = 0;
    ++ups_;
  }
  return *this;
}

template <typename bit_t>
std::pair<bit_t, bit_t> BasisNoNpIterator<bit_t>::operator*() const {
  return {ups_, dns_};
}

template <typename bit_t>
bool BasisNoNpIterator<bit_t>::operator!=(
    BasisNoNpIterator<bit_t> const &rhs) const {
  return (ups_ != rhs.ups_) || (dns_ != rhs.dns_);
}

template class BasisNoNpIterator<uint32_t>;
template class BasisNoNpIterator<uint64_t>;

} // namespace xdiag::basis::electron
