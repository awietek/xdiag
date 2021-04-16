#include "index_electron.h"

namespace hydra {

template <class bit_t, class idx_t>
IndexElectron<bit_t, idx_t>::IndexElectron(basis_t const &basis)
  : size_(basis.size()),
    n_sites_(basis.n_sites()),
    qn_(basis.qn()),
    basis_up_(n_sites_, qn_.n_up),
    basis_dn_(n_sites_, qn_.n_dn),
    index_up_(basis_up_),
    index_dn_(basis_dn_)
{ }

template class IndexElectron<uint16, int16>;
template class IndexElectron<uint16, int32>;
template class IndexElectron<uint16, int64>;
  
template class IndexElectron<uint32, int16>;
template class IndexElectron<uint32, int32>;
template class IndexElectron<uint32, int64>;
  
template class IndexElectron<uint64, int16>;
template class IndexElectron<uint64, int32>;
template class IndexElectron<uint64, int64>;
  
} // namespace hydra
