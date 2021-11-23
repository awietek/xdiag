#include "tj_indexing.h"

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/combinatorics/binomial.h>

namespace hydra::indexing {

using namespace combinatorics;

template <typename bit_t>
tJIndexing<bit_t>::tJIndexing(int n_sites, int n_up, int n_dn)
    : n_sites_(n_sites), n_up_(n_up), n_dn_(n_dn), n_spins_(n_up + n_dn),
      n_holes_(n_sites - n_spins_), size_holes_(binomial(n_sites, n_holes_)),
      size_spins_(binomial(n_spins_, n_up)), size_(size_spins_ * size_holes_),
      lintable_holes_(n_sites, n_holes_), lintable_spins_(n_spins_, n_up) {
  utils::check_nup_ndn_tj(n_sites, n_up, n_dn, "tJ");
}

template class tJIndexing<uint16_t>;
template class tJIndexing<uint32_t>;
template class tJIndexing<uint64_t>;

} // namespace hydra::indexing
