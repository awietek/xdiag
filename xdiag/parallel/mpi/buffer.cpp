// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "buffer.hpp"

namespace xdiag::mpi {

void Buffer::clean() {
  clean_send();
  clean_recv();
}

void Buffer::clean_send() { std::fill(send_.begin(), send_.end(), 0); }
void Buffer::clean_recv() { std::fill(recv_.begin(), recv_.end(), 0); }
void Buffer::swap() { std::swap(send_, recv_); }
} // namespace xdiag::mpi
