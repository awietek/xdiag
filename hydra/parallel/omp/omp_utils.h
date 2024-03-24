#pragma once

#ifdef _OPENMP

#include <hydra/extern/gsl/span>
#include <omp.h>
#include <vector>

#include <hydra/combinatorics/binomial.h>
#include <hydra/common.h>

namespace hydra::omp {
int64_t get_omp_start(int64_t size);
int64_t get_omp_end(int64_t size);
std::pair<int64_t, int64_t> get_omp_start_end(int64_t size);
std::pair<int64_t, int64_t> get_omp_subsets_start_end(int64_t n);
std::pair<int64_t, int64_t> get_omp_combinations_start_end(int64_t n,
                                                           int64_t k);

template <typename T>
inline gsl::span<T> get_omp_span(std::vector<T> const &vec) {
  auto [start, end] = get_omp_start_end(vec.size());
  return gsl::span<T>(vec.data() + start, end - start);
}

template <typename T>
inline std::vector<T>
combine_vectors(std::vector<std::vector<T>> const &vec_of_vec) {
  int64_t size = 0;
  for (auto const &vec : vec_of_vec) {
    size += vec.size();
  }

  std::vector<T> total_vec(size);
  int64_t offset = 0;

  for (auto const &vec : vec_of_vec) {

#pragma omp parallel for
    for (int64_t i = 0; i < (int64_t)vec.size(); ++i) {
      total_vec[offset + i] = vec[i];
    }
    offset += vec.size();
  }
  return total_vec;
}

template <typename T>
inline std::vector<T>
combine_vectors_copy(std::vector<std::vector<T>> const &vec_of_vec) {
  int64_t size = 0;
  for (auto const &vec : vec_of_vec) {
    size += vec.size();
  }

  std::vector<T> total_vec(size);
  int64_t offset = 0;

  for (auto const &vec : vec_of_vec) {
    std::copy(vec.begin(), vec.end(), total_vec.begin() + offset);
    offset += vec.size();
  }
  return total_vec;
}

} // namespace hydra::omp

#endif // _OPENMP
