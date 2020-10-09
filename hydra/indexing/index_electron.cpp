#include "index_electron.h"

namespace hydra {

template <class bit_t, class idx_t>
IndexElectron<bit_t, idx_t>::IndexElectron(basis_t const &basis)
  : size_(basis.size()),
    n_sites_(basis.n_sites()),
    qn_(basis.qn())
{
  BasisSpinHalf<bit_t> basis_up(n_sites_, qn_.n_up);
  BasisSpinHalf<bit_t> basis_down(n_sites_, qn_.n_dn);
  index_up_ = indexing_t(basis_up);
  index_dn_ = indexing_t(basis_down);
}

} // namespace hydra
