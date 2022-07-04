#pragma once

#include <algorithm>

#include <hydra/common.h>

namespace hydra::mpi {

class Buffer {
public:
  Buffer() = default;

  template <typename coeff_t> void reserve(idx_t size) {
    reserve<coeff_t>(size, size);
  }
  template <typename coeff_t> void reserve(idx_t size_send, idx_t size_recv);
  template <typename coeff_t> coeff_t *send();
  template <typename coeff_t> coeff_t *recv();

  void clean_send() { std::fill(send_.begin(), send_.end(), 0); }
  void clean_recv() { std::fill(recv_.begin(), recv_.end(), 0); }

private:
  std::vector<double> send_;
  std::vector<double> recv_;
};

inline Buffer buffer;

} // namespace hydra::mpi
