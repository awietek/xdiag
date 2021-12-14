#include "electron_indexing_no_np.h"

namespace hydra::indexing {

using namespace combinatorics;

template <typename bit_t>
ElectronIndexingNoNp<bit_t>::ElectronIndexingNoNp(int n_sites)
    : n_sites_(n_sites), size_ups_((idx_t)1 << n_sites),
      size_dns_((idx_t)1 << n_sites), size_(size_ups_ * size_dns_) {}

template class ElectronIndexingNoNp<uint16_t>;
template class ElectronIndexingNoNp<uint32_t>;
template class ElectronIndexingNoNp<uint64_t>;

} // namespace hydra::indexing
