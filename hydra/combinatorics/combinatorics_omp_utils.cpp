#include "combinatorics_omp_utils.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <hydra/combinatorics/binomial.h>

namespace hydra::combinatorics {

std::pair<idx_t, idx_t> get_omp_subsets_start_end(int n) {

#ifdef _OPENMP
  idx_t myid = omp_get_thread_num();
  idx_t rank = omp_get_num_threads();
  idx_t size = (idx_t)1 << n;
  idx_t chunksize = size / rank;

  idx_t start = myid * chunksize;
  idx_t end = (myid == rank - 1) ? size : (myid + 1) * chunksize;
#else
  idx_t start = 0;
  idx_t end = (idx_t)1 << n;
#endif

  return {start, end};
}

std::pair<idx_t, idx_t> get_omp_combinations_start_end(int n, int k) {

#ifdef _OPENMP
  idx_t myid = omp_get_thread_num();
  idx_t rank = omp_get_num_threads();
  idx_t size = binomial(n, k);
  idx_t chunksize = size / rank;

  idx_t start = myid * chunksize;
  idx_t end = (myid == rank - 1) ? size : (myid + 1) * chunksize;
#else
  idx_t start = 0;
  idx_t end = binomial(n, k);
#endif

  return {start, end};
}

} // namespace hydra::combinatorics
