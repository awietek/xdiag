// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "blocks.hpp"

#include <xdiag/parallel/mpi/cdot_distributed.hpp>

namespace xdiag {

int64_t dim(Block const &block) {
  return std::visit([&](auto &&b) { return b.dim(); }, block);
}

int64_t size(Block const &block) {
  return std::visit([&](auto &&b) { return b.size(); }, block);
}

int64_t nsites(Block const &block) {
  return std::visit([&](auto &&b) { return b.nsites(); }, block);
}

bool isreal(Block const &block) {
  return std::visit([&](auto &&b) { return b.isreal(); }, block);
}

std::ostream &operator<<(std::ostream &out, Block const &block) {
  std::visit([&](auto &&block) { out << block; }, block);
  return out;
}
std::string to_string(Block const &block) { return to_string_generic(block); }

} // namespace xdiag
