#ifndef HYDRA_INDEXING_INDEXTABLE_
#define HYDRA_INDEXING_INDEXTABLE_

#include <vector>

#include <hydra/common.h>

namespace hydra {

template <class basis_t, class idx_t = std_idx_t> class IndexTable {
public:
  using state_t = typename basis_t::state_t;

  IndexTable() = default;
  IndexTable(basis_t const &basis) : basis_(basis), n_sites_(basis.n_sites()) {
    indices_.resize(basis.rawsize());
    idx_t idx = 0;
    for (auto state : basis) {
      states_.push_back(state);
      indices_[rawidx(state, n_sites_)] = idx;
      ++idx;
    }
  }

  inline idx_t index(state_t const &state) const {
    return indices_[rawidx(state, n_sites_)];
  }
  inline state_t state(idx_t const &index) const { return states_[index]; }
  inline idx_t size() const { return basis_.size(); }
  inline idx_t rawsize() const { return basis_.rawsize(); }

private:
  basis_t basis_;
  number_t n_sites_;
  std::vector<idx_t> indices_;
  std::vector<state_t> states_;
};

} // namespace hydra

#endif
