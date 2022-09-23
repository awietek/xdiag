#include "indexing_np.h"

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/combinatorics/binomial.h>

namespace hydra::indexing::electron {

using namespace combinatorics;

template <typename bit_t>
IndexingNp<bit_t>::IndexingNp(int n_sites, int n_up, int n_dn)
    : n_sites_(n_sites), n_up_(n_up), n_dn_(n_dn),
      size_ups_(binomial(n_sites, n_up)), size_dns_(binomial(n_sites, n_dn)),
      size_(size_ups_ * size_dns_), lintable_ups_(n_sites, n_up),
      lintable_dns_(n_sites, n_dn) {
  utils::check_nup_ndn_electron(n_sites, n_up, n_dn, "Electron");
}

template class IndexingNp<uint16_t>;
template class IndexingNp<uint32_t>;
template class IndexingNp<uint64_t>;

} // namespace hydra::indexing::electron
