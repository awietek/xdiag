#include "blocks.hpp"

#include <xdiag/parallel/mpi/cdot_distributed.hpp>

namespace xdiag {

int64_t dim(Block const &block) {
  return std::visit([&](auto &&b) { return b.dim(); }, block);
}

int64_t size(Block const &block) {
  return std::visit([&](auto &&b) { return b.size(); }, block);
}

int64_t n_sites(Block const &block) {
  return std::visit([&](auto &&b) { return b.n_sites(); }, block);
}

bool isreal(Block const &block) {
  return std::visit([&](auto &&b) { return b.isreal(); }, block);
}

bool isdistributed(Block const &block) {
  return std::visit(
      overload{
          [&](Spinhalf const &) -> bool { return false; },
          [&](tJ const &) -> bool { return false; },
          [&](Electron const &) -> bool { return false; },
#ifdef XDIAG_USE_MPI
          [&](SpinhalfDistributed const &) -> bool { return true; },
          [&](tJDistributed const &) -> bool { return true; },
#endif
          [&](auto &&) -> bool { return false; },
      },
      block);
}

std::ostream &operator<<(std::ostream &out, Block const &block) {
  std::visit([&](auto &&block) { out << block; }, block);
  return out;
}
std::string to_string(Block const &block) { return to_string_generic(block); }

} // namespace xdiag
