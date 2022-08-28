#ifdef HYDRA_ENABLE_OPENMP

#include "omp_utils.h"

namespace hydra::omp {

idx_t get_omp_start(idx_t size) {
  idx_t myid = omp_get_thread_num();
  idx_t rank = omp_get_num_threads();
  idx_t chunksize = size / rank;
  return myid * chunksize;
}

idx_t get_omp_end(idx_t size) {
  idx_t myid = omp_get_thread_num();
  idx_t rank = omp_get_num_threads();
  idx_t chunksize = size / rank;
  return (myid == rank - 1) ? size : (myid + 1) * chunksize;
}

std::pair<idx_t, idx_t> get_omp_start_end(idx_t size) {
  idx_t myid = omp_get_thread_num();
  idx_t rank = omp_get_num_threads();
  idx_t chunksize = size / rank;

  idx_t start = myid * chunksize;
  idx_t end = (myid == rank - 1) ? size : (myid + 1) * chunksize;

  return {start, end};
}

std::pair<idx_t, idx_t> get_omp_subsets_start_end(int n) {
  idx_t size = (idx_t)1 << n;
  return get_omp_start_end(size);
}

std::pair<idx_t, idx_t> get_omp_combinations_start_end(int n, int k) {
  idx_t size = combinatorics::binomial(n, k);
  return get_omp_start_end(size);
}

} // namespace hydra::omp

#endif
