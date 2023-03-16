#include "indexing_np.h"

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/combinatorics/binomial.h>

namespace hydra::indexing::tj {

using namespace combinatorics;

template <typename bit_t>
IndexingNp<bit_t>::IndexingNp(int n_sites, int n_up, int n_dn)
    : n_sites_(n_sites), n_up_(n_up), n_dn_(n_dn),
      size_ups_(binomial(n_sites, n_up)),
      size_dncs_(binomial(n_sites - n_up, n_dn)), size_(size_ups_ * size_dncs_),
      sitesmask_(((bit_t)1 << n_sites) - 1), lintable_ups_(n_sites, n_up),
      lintable_dncs_(n_sites - n_up, n_dn) {
  utils::check_nup_ndn_tj(n_sites, n_up, n_dn, "tJ");
}

template <typename bit_t> int IndexingNp<bit_t>::n_sites() const {
  return n_sites_;
}
template <typename bit_t> int IndexingNp<bit_t>::n_up() const { return n_up_; }
template <typename bit_t> int IndexingNp<bit_t>::n_dn() const { return n_dn_; }

template <typename bit_t> idx_t IndexingNp<bit_t>::size_ups() const {
  return size_ups_;
}
template <typename bit_t> idx_t IndexingNp<bit_t>::size_dncs() const {
  return size_dncs_;
}
template <typename bit_t> idx_t IndexingNp<bit_t>::size() const {
  return size_;
}
template <typename bit_t>
idx_t IndexingNp<bit_t>::ups_offset(idx_t idx_ups) const {
  return idx_ups * size_dncs_;
}

template <typename bit_t>
Combinations<bit_t> IndexingNp<bit_t>::states_ups() const {
  return Combinations<bit_t>(n_sites_, n_up_);
}

template <typename bit_t>
Combinations<bit_t> IndexingNp<bit_t>::states_dncs(bit_t ups) const {
  (void)ups;
  assert(bitops::popcnt(ups) == n_up_);
  return Combinations<bit_t>(n_sites_ - n_up_, n_dn_);
}

template <typename bit_t>
CombinationsIndex<bit_t> IndexingNp<bit_t>::states_indices_ups() const {
  return CombinationsIndex<bit_t>(n_sites_, n_up_);
}

template <typename bit_t>
CombinationsIndex<bit_t>
IndexingNp<bit_t>::states_indices_dncs(bit_t ups) const {
  (void)ups;
  assert(bitops::popcnt(ups) == n_up_);
  return CombinationsIndex<bit_t>(n_sites_ - n_up_, n_dn_);
}

#ifdef _OPENMP
template <typename bit_t>
CombinationsThread<bit_t> IndexingNp<bit_t>::states_ups_thread() const {
  return CombinationsThread<bit_t>(n_sites_, n_up_);
}

template <typename bit_t>
CombinationsThread<bit_t>
IndexingNp<bit_t>::states_dncs_thread(bit_t ups) const {
  (void)ups;
  assert(bitops::popcnt(ups) == n_up_);
  return CombinationsThread<bit_t>(n_sites_ - n_up_, n_dn_);
}

template <typename bit_t>
CombinationsIndexThread<bit_t>
IndexingNp<bit_t>::states_indices_ups_thread() const {
  return CombinationsIndexThread<bit_t>(n_sites_, n_up_);
}

template <typename bit_t>
CombinationsIndexThread<bit_t>
IndexingNp<bit_t>::states_indices_dncs_thread(bit_t ups) const {
  (void)ups;
  assert(bitops::popcnt(ups) == n_up_);
  return CombinationsIndexThread<bit_t>(n_sites_ - n_up_, n_dn_);
}
#endif

template class IndexingNp<uint16_t>;
template class IndexingNp<uint32_t>;
template class IndexingNp<uint64_t>;

} // namespace hydra::indexing::tj
