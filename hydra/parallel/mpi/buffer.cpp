#include "buffer.h"

namespace hydra::mpi {

void Buffer::clean() {
  clean_send();
  clean_recv();
}

void Buffer::clean_send() { std::fill(send_.begin(), send_.end(), 0); }
void Buffer::clean_recv() { std::fill(recv_.begin(), recv_.end(), 0); }
void Buffer::swap() { std::swap(send_, recv_); }
} // namespace hydra::mpi
