#pragma once

#include <vector>

#include <hydra/common.h>
#include <hydra/utils/bitops.h>

namespace hydra {

template <class bit_t = std_bit_t> class LinTable {
public:
  LinTable() = default;
  LinTable(int n, int k);

  inline idx_t index(bit_t bits) const {
    return left_indices_[bits >> n_right_] +
           right_indices_[utils::gbits(bits, n_right_, 0)];
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

} // namespace hydra
