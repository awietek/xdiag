#include "basis_np.hpp"

#include <xdiag/combinatorics/binomial.hpp>

namespace xdiag::basis::tj {

using namespace combinatorics;

template <typename bit_t>
BasisNp<bit_t>::BasisNp(int64_t nsites, int64_t nup, int64_t ndn) try
    : nsites_(nsites), nup_(nup), ndn_(ndn),
      size_ups_(binomial(nsites, nup)),
      size_dncs_(binomial(nsites - nup, ndn)), size_(size_ups_ * size_dncs_),
      sitesmask_(((bit_t)1 << nsites) - 1), lintable_ups_(nsites, nup),
      lintable_dncs_(nsites - nup, ndn) {
  check_nsites_work_with_bits<bit_t>(nsites_);

  if (nsites < 0) {
    XDIAG_THROW("nsites < 0");
  } else if ((nup < 0) || (ndn < 0)) {
    XDIAG_THROW("nup < 0 or ndn < 0");
  } else if ((nup + ndn) > nsites) {
    XDIAG_THROW("nup + ndn > nsites");
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename bit_t> int64_t BasisNp<bit_t>::nsites() const {
  return nsites_;
}
template <typename bit_t> int64_t BasisNp<bit_t>::nup() const { return nup_; }
template <typename bit_t> int64_t BasisNp<bit_t>::ndn() const { return ndn_; }

template <typename bit_t> int64_t BasisNp<bit_t>::size_ups() const {
  return size_ups_;
}
template <typename bit_t> int64_t BasisNp<bit_t>::size_dncs() const {
  return size_dncs_;
}
template <typename bit_t> int64_t BasisNp<bit_t>::size() const { return size_; }
template <typename bit_t> int64_t BasisNp<bit_t>::dim() const { return size_; }
template <typename bit_t>
int64_t BasisNp<bit_t>::ups_offset(int64_t idx_ups) const {
  return idx_ups * size_dncs_;
}
template <typename bit_t>
typename BasisNp<bit_t>::iterator_t BasisNp<bit_t>::begin() const {
  return iterator_t(nsites_, nup_, ndn_, true);
}
template <typename bit_t>
typename BasisNp<bit_t>::iterator_t BasisNp<bit_t>::end() const {
  return iterator_t(nsites_, nup_, ndn_, false);
}

template <typename bit_t>
Combinations<bit_t> BasisNp<bit_t>::states_ups() const {
  return Combinations<bit_t>(nsites_, nup_);
}

template <typename bit_t>
Combinations<bit_t> BasisNp<bit_t>::states_dncs(bit_t ups) const {
  (void)ups;
  assert(bits::popcnt(ups) == nup_);
  return Combinations<bit_t>(nsites_ - nup_, ndn_);
}

template <typename bit_t>
CombinationsIndex<bit_t> BasisNp<bit_t>::states_indices_ups() const {
  return CombinationsIndex<bit_t>(nsites_, nup_);
}

template <typename bit_t>
CombinationsIndex<bit_t> BasisNp<bit_t>::states_indices_dncs(bit_t ups) const {
  (void)ups;
  assert(bits::popcnt(ups) == nup_);
  return CombinationsIndex<bit_t>(nsites_ - nup_, ndn_);
}

#ifdef _OPENMP
template <typename bit_t>
CombinationsThread<bit_t> BasisNp<bit_t>::states_ups_thread() const {
  return CombinationsThread<bit_t>(nsites_, nup_);
}

template <typename bit_t>
CombinationsThread<bit_t> BasisNp<bit_t>::states_dncs_thread(bit_t ups) const {
  (void)ups;
  assert(bits::popcnt(ups) == nup_);
  return CombinationsThread<bit_t>(nsites_ - nup_, ndn_);
}

template <typename bit_t>
CombinationsIndexThread<bit_t>
BasisNp<bit_t>::states_indices_ups_thread() const {
  return CombinationsIndexThread<bit_t>(nsites_, nup_);
}

template <typename bit_t>
CombinationsIndexThread<bit_t>
BasisNp<bit_t>::states_indices_dncs_thread(bit_t ups) const {
  (void)ups;
  assert(bits::popcnt(ups) == nup_);
  return CombinationsIndexThread<bit_t>(nsites_ - nup_, ndn_);
}
#endif

template class BasisNp<uint32_t>;
template class BasisNp<uint64_t>;

template <typename bit_t>
BasisNpIterator<bit_t>::BasisNpIterator(int64_t nsites, int64_t nup,
                                        int64_t ndn, bool begin)
    : sitesmask_(((bit_t)1 << nsites) - 1) {
  begindns_ = ((bit_t)1 << ndn) - 1;
  bit_t beginups = ((bit_t)1 << nup) - 1;
  bit_t end_ups = get_next_pattern(beginups << (nsites - nup));
  end_dns_ = get_next_pattern(begindns_ << (nsites - nup - ndn));
  ups_ = begin ? beginups : end_ups;
  dns_ = begindns_;
}

template <typename bit_t>
BasisNpIterator<bit_t> &BasisNpIterator<bit_t>::operator++() {
  dns_ = get_next_pattern(dns_);
  if (dns_ == end_dns_) {
    dns_ = begindns_;
    ups_ = get_next_pattern(ups_);
  }
  return *this;
}

template <typename bit_t>
std::pair<bit_t, bit_t> BasisNpIterator<bit_t>::operator*() const {
  bit_t not_ups = (~ups_) & sitesmask_;
  return {ups_, bits::deposit(dns_, not_ups)};
}

template <typename bit_t>
bool BasisNpIterator<bit_t>::operator!=(
    BasisNpIterator<bit_t> const &rhs) const {
  return (ups_ != rhs.ups_) || (dns_ != rhs.dns_);
}

template class BasisNpIterator<uint32_t>;
template class BasisNpIterator<uint64_t>;

} // namespace xdiag::basis::tj
