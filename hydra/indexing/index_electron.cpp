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

template class IndexElectron<uint16, uint16>;
template class IndexElectron<uint16, uint32>;
template class IndexElectron<uint16, uint64>;
  
template class IndexElectron<uint32, uint16>;
template class IndexElectron<uint32, uint32>;
template class IndexElectron<uint32, uint64>;
  
template class IndexElectron<uint64, uint16>;
template class IndexElectron<uint64, uint32>;
template class IndexElectron<uint64, uint64>;
  
} // namespace hydra
