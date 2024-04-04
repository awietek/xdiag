#pragma once
#ifdef XDIAG_USE_MPI

#include <xdiag/common.h>
#include <mpi.h>

#include <algorithm>

namespace xdiag::mpi {

class Buffer {
public:
  Buffer() = default;

  template <typename coeff_t> inline void reserve(int64_t size) {
    reserve<coeff_t>(size, size);
  }

  template <typename coeff_t>
  inline void reserve(int64_t size_send, int64_t size_recv) {
    int coeff_size = sizeof(coeff_t);
    int char_size = sizeof(char);
    int fac = (coeff_size > char_size) ? coeff_size / char_size : 1;

    if (fac * size_send > (int64_t)send_.size()) {
      send_.resize(fac * size_send);
    }
    if (fac * size_recv > (int64_t)recv_.size()) {
      recv_.resize(fac * size_recv);
    }
  }
  template <typename coeff_t> inline coeff_t *send() {
    return reinterpret_cast<coeff_t *>(send_.data());
  }
  template <typename coeff_t> inline coeff_t *recv() {
    return reinterpret_cast<coeff_t *>(recv_.data());
  }

  void clean();
  void clean_send();
  void clean_recv();
  void swap();
private:
  std::vector<char> send_;
  std::vector<char> recv_;
};

inline Buffer buffer;

} // namespace xdiag::mpi
#endif
