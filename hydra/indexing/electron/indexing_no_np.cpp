#include "indexing_no_np.h"

namespace hydra::indexing::electron {

using namespace combinatorics;

template <typename bit_t>
IndexingNoNp<bit_t>::IndexingNoNp(int n_sites)
    : n_sites_(n_sites), size_ups_((idx_t)1 << n_sites),
      size_dns_((idx_t)1 << n_sites), size_(size_ups_ * size_dns_) {}

template class IndexingNoNp<uint16_t>;
template class IndexingNoNp<uint32_t>;
template class IndexingNoNp<uint64_t>;

} // namespace hydra::indexing::electron
