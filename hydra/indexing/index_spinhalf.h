#ifndef HYDRA_INDEXING_INDEXSPINHALF_
#define HYDRA_INDEXING_INDEXSPINHALF_

#include <hydra/bases/basis_spinhalf.h>
#include <hydra/common.h>
#include <hydra/states/state_spinhalf.h>

namespace hydra {

template <class bit_t = std_bit_t, class idx_t = std_idx_t>
class IndexSpinHalf {
public:
  using state_t = state_spinhalf<bit_t>;
  using qn_t = qn_spinhalf;
  using basis_t = BasisSpinHalf<bit_t>;

  IndexSpinHalf() = default;
  IndexSpinHalf(const basis_t &basis);

  idx_t index(const state_t &state) const;
  state_t state(const idx_t &index) const;
  idx_t size() const { return size_; }

private:
  idx_t size_;
  number_t n_sites_;
  qn_t qn_;
};

} // namespace hydra
#endif
