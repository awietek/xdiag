#ifndef HYDRA_INDEXING_INDEXTABLE_
#define HYDRA_INDEXING_INDEXTABLE_

#include <vector>

#include <hydra/common.h>

namespace hydra {

template <class hilbertspace_t, class idx_t = std_idx_t> class IndexTable {
public:
  using state_t = typename hilbertspace_t::state_t;

  IndexTable() = default;
  IndexTable(hilbertspace_t const &hilbertspace);

  inline idx_t index(state_t const &state) const {
    return indices_[rawidx(state, n_sites_)];
  }
  inline state_t state(idx_t const &index) const { return states_[index]; }
  inline idx_t size() const { return size_; }

private:
  idx_t size_;
  number_t n_sites_;
  std::vector<idx_t> indices_;
  std::vector<state_t> states_;
};

} // namespace hydra

#endif
