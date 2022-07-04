#include "buffer.h"

namespace hydra::mpi {

template <> void Buffer::reserve<double>(idx_t size_send, idx_t size_recv) {
  if (size_send > send_.size()) {
    send_.resize(size_send);
  }
  if (size_recv > recv_.size()) {
    recv_.resize(size_recv);
  }
}

template <> void Buffer::reserve<complex>(idx_t size_send, idx_t size_recv) {
  if (2 * size_send > send_.size()) {
    send_.resize(2 * size_send);
  }
  if (2 * size_recv > recv_.size()) {
    recv_.resize(2 * size_recv);
  }
}

template <> double *Buffer::send<double>() { return send_.data(); }
template <> complex *Buffer::send<complex>() {
  return reinterpret_cast<complex *>(send_.data());
}

template <> double *Buffer::recv<double>() { return recv_.data(); }
template <> complex *Buffer::recv<complex>() {
  return reinterpret_cast<complex *>(recv_.data());
}

} // namespace hydra::mpi
