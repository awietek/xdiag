#ifdef _OPENMP

#include "omp_utils.h"

namespace hydra::omp {

int64_t get_omp_start(int64_t size) {
  int64_t myid = omp_get_thread_num();
  int64_t rank = omp_get_num_threads();
  int64_t chunksize = size / rank;
  return myid * chunksize;
}

int64_t get_omp_end(int64_t size) {
  int64_t myid = omp_get_thread_num();
  int64_t rank = omp_get_num_threads();
  int64_t chunksize = size / rank;
  return (myid == rank - 1) ? size : (myid + 1) * chunksize;
}

std::pair<int64_t, int64_t> get_omp_start_end(int64_t size) {
  int64_t myid = omp_get_thread_num();
  int64_t rank = omp_get_num_threads();
  int64_t chunksize = size / rank;

  int64_t start = myid * chunksize;
  int64_t end = (myid == rank - 1) ? size : (myid + 1) * chunksize;

  return {start, end};
}

std::pair<int64_t, int64_t> get_omp_subsets_start_end(int64_t n) {
  int64_t size = (int64_t)1 << n;
  return get_omp_start_end(size);
}

std::pair<int64_t, int64_t> get_omp_combinations_start_end(int64_t n, int64_t k) {
  int64_t size = combinatorics::binomial(n, k);
  return get_omp_start_end(size);
}

} // namespace hydra::omp

#endif
