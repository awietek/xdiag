#include "spinhalf_indexing_sz.h"

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/combinatorics/binomial.h>

namespace hydra::indexing {

template <typename bit_t>
SpinhalfIndexingSz<bit_t>::SpinhalfIndexingSz(int n_sites, int nup)
    : n_sites_(n_sites), n_up_(nup), lintable_(n_sites, nup),
      size_(combinatorics::binomial(n_sites, nup)) {
  utils::check_nup_spinhalf(n_sites, nup, "Spinhalf");
}

template class SpinhalfIndexingSz<uint16_t>;
template class SpinhalfIndexingSz<uint32_t>;
template class SpinhalfIndexingSz<uint64_t>;

} // namespace hydra::indexing
