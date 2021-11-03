#pragma once

#include <vector>

#include <hydra/common.h>
#include <hydra/bitops/bitops.h>

namespace hydra::indexing {

template <class bit_t = std_bit_t> class LinTable {
public:
  LinTable() = default;
  LinTable(int n, int k);

  inline idx_t index(bit_t bits) const {
    return left_indices_[bits >> n_right_] +
           right_indices_[bitops::gbits(bits, 0, n_right_)];
  }

private:
  int n_, k_;
  int n_left_;
  int n_right_;

  idx_t left_table_size_;
  idx_t right_table_size_;

  std::vector<idx_t> left_indices_;
  std::vector<idx_t> right_indices_;
};

} // namespace hydra::indexing
