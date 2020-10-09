#ifndef HYDRA_INDEXING_LINTABLE_
#define HYDRA_INDEXING_LINTABLE_

#include <vector>

#include <hydra/common.h>
#include <hydra/utils/bitops.h>

namespace hydra {

using utils::gbits;

template <class hilbertspace_t, class idx_t = std_idx_t> class LinTable {
public:
  using state_t = typename hilbertspace_t::state_t;

  LinTable() = default;
  LinTable(hilbertspace_t const &hilbertspace);

  idx_t index(state_t const &state) const;
  inline state_t state(idx_t const &index) const { return states_[index]; }

  inline idx_t size() const { return total_size_; }

private:
  idx_t total_size_;
  idx_t right_size_;
  idx_t left_size_;
  int site_divider_;
  std::vector<idx_t> right_indices_;
  std::vector<idx_t> left_indices_;
  int num_of_right_sites;
  int num_of_left_sites;
  std::vector<state_t> states_;
};

} // namespace hydra

#endif
