#pragma once

#ifdef HYDRA_ENABLE_OPENMP

#include <hydra/combinatorics/binomial.h>
#include <hydra/common.h>
#include <lila/external/gsl/span>
#include <omp.h>
#include <vector>

namespace hydra::omp {
idx_t get_omp_start(idx_t size);
idx_t get_omp_end(idx_t size);
std::pair<idx_t, idx_t> get_omp_start_end(idx_t size);
std::pair<idx_t, idx_t> get_omp_subsets_start_end(int n);
std::pair<idx_t, idx_t> get_omp_combinations_start_end(int n, int k);

template <typename T>
inline gsl::span<T> get_omp_span(std::vector<T> const &vec) {
  auto [start, end] = get_omp_start_end(vec.size());
  return gsl::span<T>(vec.data() + start, end - start);
}

template <typename T>
inline std::vector<T>
combine_vectors(std::vector<std::vector<T>> const &vec_of_vec) {
  idx_t size = 0;
  for (auto const &vec : vec_of_vec) {
    size += vec.size();
  }

  std::vector<T> total_vec(size);
  idx_t offset = 0;

  for (auto const &vec : vec_of_vec) {

#pragma openmp parallel for
    for (idx_t i = 0; i < (idx_t)vec.size(); ++i) {
      total_vec[offset + i] = vec[i];
    }
    offset += vec.size();
  }
  return total_vec;
}

template <typename T>
inline std::vector<T>
combine_vectors_copy(std::vector<std::vector<T>> const &vec_of_vec) {
  idx_t size = 0;
  for (auto const &vec : vec_of_vec) {
    size += vec.size();
  }

  std::vector<T> total_vec(size);
  idx_t offset = 0;

  for (auto const &vec : vec_of_vec) {
    std::copy(vec.begin(), vec.end(), total_vec.begin() + offset);
    offset += vec.size();
  }
  return total_vec;
}

} // namespace hydra::omp

#endif // HYDRA_ENABLE_OPENMP
