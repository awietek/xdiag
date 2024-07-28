#include "basis_np.hpp"
#include <xdiag/combinatorics/binomial.hpp>
#include <xdiag/combinatorics/bit_patterns.hpp>

namespace xdiag::basis::electron {

using namespace combinatorics;

template <typename bit_t>
BasisNp<bit_t>::BasisNp(int n_sites, int n_up, int n_dn)
    : n_sites_(n_sites), n_up_(n_up), n_dn_(n_dn),
      size_ups_(binomial(n_sites, n_up)), size_dns_(binomial(n_sites, n_dn)),
      size_(size_ups_ * size_dns_), lintable_ups_(n_sites, n_up),
      lintable_dns_(n_sites, n_dn) {
  if (n_sites < 0) {
    throw(std::invalid_argument("n_sites < 0"));
  } else if ((n_up < 0) || (n_up > n_sites)) {
    throw(std::invalid_argument("Invalid value of nup"));
  } else if ((n_dn < 0) || (n_dn > n_sites)) {
    throw(std::invalid_argument("Invalid value of ndn"));
  }
}

template <typename bit_t> int BasisNp<bit_t>::n_sites() const {
  return n_sites_;
}
template <typename bit_t> int BasisNp<bit_t>::n_up() const { return n_up_; }
template <typename bit_t> int BasisNp<bit_t>::n_dn() const { return n_dn_; }

template <typename bit_t> int64_t BasisNp<bit_t>::size_ups() const {
  return size_ups_;
}
template <typename bit_t> int64_t BasisNp<bit_t>::size_dns() const {
  return size_dns_;
}

template <typename bit_t> int64_t BasisNp<bit_t>::dim() const { return size_; }
template <typename bit_t> int64_t BasisNp<bit_t>::size() const { return size_; }
template <typename bit_t>
typename BasisNp<bit_t>::iterator_t BasisNp<bit_t>::begin() const {
  return iterator_t(n_sites_, n_up_, n_dn_, true);
}
template <typename bit_t>
typename BasisNp<bit_t>::iterator_t BasisNp<bit_t>::end() const {
  return iterator_t(n_sites_, n_up_, n_dn_, false);
}

template <typename bit_t>
Combinations<bit_t> BasisNp<bit_t>::states_ups() const {
  return Combinations<bit_t>(n_sites_, n_up_);
}

template <typename bit_t>
Combinations<bit_t> BasisNp<bit_t>::states_dns() const {
  return Combinations<bit_t>(n_sites_, n_dn_);
}

template <typename bit_t>
CombinationsIndex<bit_t> BasisNp<bit_t>::states_indices_ups() const {
  return CombinationsIndex<bit_t>(n_sites_, n_up_);
}

template <typename bit_t>
CombinationsIndex<bit_t> BasisNp<bit_t>::states_indices_dns() const {
  return CombinationsIndex<bit_t>(n_sites_, n_dn_);
}

#ifdef _OPENMP
template <typename bit_t>
CombinationsThread<bit_t> BasisNp<bit_t>::states_ups_thread() const {
  return CombinationsThread<bit_t>(n_sites_, n_up_);
}

template <typename bit_t>
CombinationsThread<bit_t> BasisNp<bit_t>::states_dns_thread() const {
  return CombinationsThread<bit_t>(n_sites_, n_dn_);
}

template <typename bit_t>
CombinationsIndexThread<bit_t>
BasisNp<bit_t>::states_indices_ups_thread() const {
  return CombinationsIndexThread<bit_t>(n_sites_, n_up_);
}

template <typename bit_t>
CombinationsIndexThread<bit_t>
BasisNp<bit_t>::states_indices_dns_thread() const {
  return CombinationsIndexThread<bit_t>(n_sites_, n_dn_);
}
#endif

template class BasisNp<uint32_t>;
template class BasisNp<uint64_t>;

template <typename bit_t>
BasisNpIterator<bit_t>::BasisNpIterator(int64_t n_sites, int64_t n_up,
                                        int64_t n_dn, bool begin) {
  begin_dns_ = ((bit_t)1 << n_dn) - 1;
  bit_t begin_ups = ((bit_t)1 << n_up) - 1;
  bit_t end_ups = get_next_pattern(begin_ups << (n_sites - n_up));
  end_dns_ = get_next_pattern(begin_dns_ << (n_sites - n_dn));
  ups_ = begin ? begin_ups : end_ups;
  dns_ = begin_dns_;
}

template <typename bit_t>
BasisNpIterator<bit_t> &BasisNpIterator<bit_t>::operator++() {
  dns_ = get_next_pattern(dns_);
  if (dns_ == end_dns_) {
    dns_ = begin_dns_;
    ups_ = get_next_pattern(ups_);
  }
  return *this;
}

template <typename bit_t>
std::pair<bit_t, bit_t> BasisNpIterator<bit_t>::operator*() const {
  return {ups_, dns_};
}

template <typename bit_t>
bool BasisNpIterator<bit_t>::operator!=(
    BasisNpIterator<bit_t> const &rhs) const {
  return (ups_ != rhs.ups_) || (dns_ != rhs.dns_);
}

template class BasisNpIterator<uint32_t>;
template class BasisNpIterator<uint64_t>;

} // namespace xdiag::basis::electron
