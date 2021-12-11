#include "spinhalf.h"

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/combinatorics/binomial.h>

namespace hydra {

template <typename bit_t>
Spinhalf<bit_t>::Spinhalf(int n_sites)
    : n_sites_(n_sites), sz_conserved_(false), n_up_(invalid_n),
      n_dn_(invalid_n), sz_(invalid_n), size_(pow(2, n_sites)),
      indexing_sz_conserved_(),
      indexing_sz_not_conserved_(std::make_shared<idxng_no_sz_t>(n_sites)) {
  assert(n_sites >= 0);
}

template <typename bit_t>
Spinhalf<bit_t>::Spinhalf(int n_sites, int n_up)
    : n_sites_(n_sites), sz_conserved_(true), n_up_(n_up),
      n_dn_(n_sites - n_up), sz_(n_up_ - n_dn_),
      size_(combinatorics::binomial(n_sites, n_up)),
      indexing_sz_conserved_(std::make_shared<idxng_sz_t>(n_sites, n_up)),
      indexing_sz_not_conserved_() {
  utils::check_nup_spinhalf(n_sites, n_up, "Spinhalf");
}

template <typename bit_t>
bool Spinhalf<bit_t>::operator==(Spinhalf<bit_t> const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (sz_conserved_ == rhs.sz_conserved_) &&
         (sz_ == rhs.sz_) && (n_up_ == rhs.n_up_) && (n_dn_ == rhs.n_dn_);
}

template <typename bit_t>
bool Spinhalf<bit_t>::operator!=(Spinhalf<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template <typename bit_t>
indexing::SpinhalfIndexingNoSz<bit_t> const &
Spinhalf<bit_t>::indexing_sz_not_conserved() const {
  if (sz_conserved_)
    lila::Log.err("Error: wrong indexing required");
  return *indexing_sz_not_conserved_;
}

template <typename bit_t>
indexing::SpinhalfIndexing<bit_t> const &
Spinhalf<bit_t>::indexing_sz_conserved() const {
  if (!sz_conserved_)
    lila::Log.err("Error: wrong indexing required");
  return *indexing_sz_conserved_;
}

template class Spinhalf<uint16_t>;
template class Spinhalf<uint32_t>;
template class Spinhalf<uint64_t>;

} // namespace hydra
