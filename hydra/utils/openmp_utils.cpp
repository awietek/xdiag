#include "openmp_utils.h"
#ifdef _OPENMP
#include <omp.h>
#endif

namespace hydra::utils {

std::pair<idx_t, idx_t> get_omp_start_end(int size) {
#ifdef _OPENMP
  idx_t myid = omp_get_thread_num();
  idx_t rank = omp_get_num_threads();
  idx_t chunksize = size / rank;

  idx_t start = myid * chunksize;
  idx_t end = (myid == rank - 1) ? size : (myid + 1) * chunksize;
#else
  idx_t start = 0;
  idx_t end = size;
#endif
  return {start, end};
}
  
}  // namespace hydra::utils
