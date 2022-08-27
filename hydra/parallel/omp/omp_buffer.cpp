#ifdef HYDRA_ENABLE_OPENMP
#include "omp_buffer.h"

namespace hydra::omp {

template <typename coeff_t>
OmpBuffer<coeff_t>::OmpBuffer()
    : myid_(omp_get_thread_num()), rank_(omp_get_num_threads()), start_(0),
      end_(0), size_(0), sorted_(false) {}

template <typename coeff_t> void OmpBuffer<coeff_t>::reserve(idx_t size) {
  std::tie(start_, end_) = get_omp_start_end(size);
  size_ = end_ - start_;
  in_out_coeffs.resize(size);
}

template class OmpBuffer<double>;
template class OmpBuffer<complex>;

} // namespace hydra::omp

#endif
