// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#ifdef _OPENMP

#include <vector>

#include <omp.h>

#include <xdiag/combinatorics/binomial.hpp>
#include <xdiag/common.hpp>
#include <xdiag/extern/gsl/span>

namespace xdiag::omp {
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
} // namespace xdiag::omp

#endif // _OPENMP

namespace xdiag::omp {
  
// The structs and definitions below are used to extract the
// number of parameters of a lambda function
// This is a bit black magic, but described here:
// https://stackoverflow.com/questions/54389831/count-the-number-of-arguments-in-a-lambda

struct any_argument {
  template <typename T> operator T &&() const;
};

template <typename Lambda, typename Is, typename = void>
struct can_accept_impl : std::false_type {};

template <typename Lambda, std::size_t... Is>
struct can_accept_impl<
    Lambda, std::index_sequence<Is...>,
    decltype(std::declval<Lambda>()(((void)Is, any_argument{})...), void())>
    : std::true_type {};

template <typename Lambda, std::size_t N>
struct can_accept : can_accept_impl<Lambda, std::make_index_sequence<N>> {};

template <typename Lambda, std::size_t Max, std::size_t N, typename = void>
struct lambda_details_impl : lambda_details_impl<Lambda, Max, N - 1> {};

template <typename Lambda, std::size_t Max, std::size_t N>
struct lambda_details_impl<Lambda, Max, N,
                           std::enable_if_t<can_accept<Lambda, N>::value>> {
  static constexpr bool is_variadic = (N == Max);
  static constexpr std::size_t argument_count = N;
};

template <typename Lambda, std::size_t Max = 50>
struct lambda_details : lambda_details_impl<Lambda, Max, Max> {};

} // namespace xdiag::omp

// This macro optionally adds the thread number to the fill call
#ifdef _OPENMP
#define XDIAG_FILL(idx_in, idx_out, coeff)                                     \
  if constexpr (omp::lambda_details<decltype(fill)>::argument_count == 3) {    \
    fill(idx_in, idx_out, coeff);                                              \
  } else {                                                                     \
    fill(idx_in, idx_out, coeff, num_thread);                                  \
  }
#else
#define XDIAG_FILL(idx_in, idx_out, coeff) fill(idx_in, idx_out, coeff)
#endif
