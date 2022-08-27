#pragma once
#ifdef HYDRA_ENABLE_OPENMP

#include <omp.h>
#include <hydra/common.h>
#include <hydra/parallel/omp/omp_utils.h>

namespace hydra::omp {

template <typename coeff_t> class OmpBuffer {
public:
  using buffer_t = std::vector<std::tuple<idx_t, idx_t, coeff_t>>;
  using iterator_t = buffer_t::iterator;
  
  OmpBuffer();

  void reserve(idx_t size);
  void append(idx_t idx_in, idx_t idx_out, coeff_t coeff);
  void sort();

  iterator_t begin() const;
  iterator_t end() const;

  void clean();

private:
  int myid_;
  int rank_;
  idx_t start_;
  idx_t end_;
  idx_t size_;
  bool sorted_;
  
  buffer_t in_out_coeffs_;
};

} // namespace hydra::omp

#endif
