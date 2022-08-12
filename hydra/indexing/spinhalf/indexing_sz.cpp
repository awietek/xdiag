#include "indexing_sz.h"

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/combinatorics/binomial.h>

namespace hydra::indexing::spinhalf {

template <typename bit_t>
IndexingSz<bit_t>::IndexingSz(int n_sites, int nup)
    : n_sites_(n_sites), n_up_(nup), lintable_(n_sites, nup),
      size_(combinatorics::binomial(n_sites, nup)), begin_(n_sites, nup, 0),
      end_(n_sites, nup, size_) {
  utils::check_nup_spinhalf(n_sites, nup, "Spinhalf");
}

template class IndexingSz<uint16_t>;
template class IndexingSz<uint32_t>;
template class IndexingSz<uint64_t>;

} // namespace hydra::indexing::spinhalf
