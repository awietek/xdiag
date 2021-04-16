#ifndef HYDRA_INDEXING_INDEXSEARCH_
#define HYDRA_INDEXING_INDEXSEARCH_

#include <vector>

#include <hydra/common.h>

namespace hydra {

template <class basis_t, class idx_t = std_idx_t> class IndexSearch {
public:
  using state_t = typename basis_t::state_t;

  IndexSearch() = default;
  IndexSearch(basis_t const &basis) : basis_(basis) {
    for (auto state : basis)
      states_.push_back(state);
  }

  idx_t index(state_t const &state) const {
    // Do binary search
    return std::lower_bound(states_.begin(), states_.end(), state) -
           states_.begin();
  }  
  state_t state(idx_t const &index) const { return states_[index]; }
  idx_t size() const { return basis_.size(); }

private:
  basis_t basis_;
  std::vector<state_t> states_;
};

} // namespace hydra

#endif
