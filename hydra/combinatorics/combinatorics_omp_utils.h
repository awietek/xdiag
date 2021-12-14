#pragma once

#include <hydra/common.h>
#include <utility>

namespace hydra::combinatorics {

std::pair<idx_t, idx_t> get_omp_subsets_start_end(int n);
std::pair<idx_t, idx_t> get_omp_combinations_start_end(int n, int k);

} // namespace hydra::combinatorics
