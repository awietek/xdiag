#ifndef HYDRA_INDEXING_INDEXSEARCH_
#define HYDRA_INDEXING_INDEXSEARCH_

#include <vector>

#include <hydra/common.h>

namespace hydra {

template <class basis_t, class idx_t = std_idx_t> class IndexSearch {
public:
  using state_t = typename basis_t::state_t;

  IndexSearch() = default;
  IndexSearch(const basis_t &basis);

  idx_t index(const state_t &state) const;
  state_t state(const idx_t &index) const { return states_[index]; }
  idx_t size() const { return size_; }

private:
  idx_t size_;
  std::vector<state_t> states_;
};

} // namespace hydra

#endif
