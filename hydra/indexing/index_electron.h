#ifndef HYDRA_INDEXING_INDEXELECTRON_
#define HYDRA_INDEXING_INDEXELECTRON_

#include <vector>

#include <hydra/bases/basis_electron.h>
#include <hydra/common.h>
#include <hydra/indexing/index_spinhalf.h>

namespace hydra {

template <class bit_t = std_bit_t, class idx_t = std_idx_t>
class IndexElectron {
public:
  using state_t = state_electron<bit_t>;
  using basis_t = BasisElectron<bit_t>;
  using qn_t = qn_electron;

  IndexElectron() = default;
  IndexElectron(const basis_t &basis);

  inline idx_t index(const state_t &state) const {
    return index_up_.index({state.ups}) * index_dn_.size() +
      index_dn_.index({state.dns});
  }

  inline state_t state(const idx_t &index) const {
    idx_t up_idx = index / index_dn_.size();
    idx_t dn_idx = index % index_dn_.size();
    return {index_up_.state(up_idx).spins, index_dn_.state(dn_idx).spins};
  }
  idx_t size() const { return size_; }

private:
  idx_t size_;
  number_t n_sites_;
  qn_t qn_;
  IndexSpinHalf<bit_t, idx_t> index_up_, index_dn_;
};

} // namespace hydra

#endif
