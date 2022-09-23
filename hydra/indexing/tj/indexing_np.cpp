#include "indexing_np.h"

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/combinatorics/binomial.h>

namespace hydra::indexing::tj {

using namespace combinatorics;

template <typename bit_t>
IndexingNp<bit_t>::IndexingNp(int n_sites, int n_up, int n_dn)
    : n_sites_(n_sites), n_up_(n_up), n_dn_(n_dn), n_spins_(n_up + n_dn),
      n_holes_(n_sites - n_spins_), size_holes_(binomial(n_sites, n_holes_)),
      size_spins_(binomial(n_spins_, n_up)), size_(size_spins_ * size_holes_),
      lintable_holes_(n_sites, n_holes_), lintable_spins_(n_spins_, n_up) {
  utils::check_nup_ndn_tj(n_sites, n_up, n_dn, "tJ");
}

template class IndexingNp<uint16_t>;
template class IndexingNp<uint32_t>;
template class IndexingNp<uint64_t>;

} // namespace hydra::indexing::tj
