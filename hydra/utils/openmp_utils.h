#pragma once

// #include <execution>

#include <hydra/common.h>

namespace hydra::utils {
template <typename T>
std::vector<T> combine_vectors(std::vector<std::vector<T>> const &vec_of_vec) {
  idx_t size = 0;
  for (auto const &vec : vec_of_vec) {
    size += vec.size();
  }

  std::vector<T> total_vec(size);
  idx_t offset = 0;

  for (auto const &vec : vec_of_vec) {
#pragma omp parallel for
    for (idx_t i = 0; i < vec.size(); ++i) {
      total_vec[offset + i] = vec[i];
    }
    offset += vec.size();
  }
  return total_vec;
}

template <typename T>
std::vector<T>
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

} // namespace hydra::utils
