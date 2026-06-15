// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

namespace xdiag::utils {

// The contiguous chunk [begin, end) of a container assigned to one OpenMP
// thread, together with the global linear index of `begin`.
template <typename iterator_t> struct ThreadRange {
  iterator_t begin;
  iterator_t end;
  int64_t index;
};

// Split a container of `container.size()` elements into `nthreads` near-equal
// contiguous chunks and return the chunk owned by thread `num_thread` (the last
// thread takes the remainder), with the global index of its first element. The
// container must provide size() and a random-access begin()/end(). With
// num_thread = 0, nthreads = 1 it returns the whole container (the serial case).
// Use as:
//
//   #pragma omp parallel
//   {
//     int num_thread = omp_get_thread_num();
//     auto rng = thread_range(container, num_thread, omp_get_num_threads());
//     int64_t idx = rng.index;
//     for (auto it = rng.begin; it != rng.end; ++it, ++idx) { ... }
//   }
template <typename container_t>
inline auto thread_range(container_t const &container, int num_thread,
                         int nthreads)
    -> ThreadRange<decltype(container.begin())> {
  int64_t size = container.size();
  int64_t index = num_thread * (size / nthreads);
  return {container.begin() + index,
          (num_thread == nthreads - 1)
              ? container.end()
              : container.begin() + (num_thread + 1) * (size / nthreads),
          index};
}

} // namespace xdiag::utils
