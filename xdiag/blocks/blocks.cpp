#include "blocks.hpp"

#include <xdiag/parallel/mpi/cdot_distributed.hpp>

namespace xdiag {

int64_t dim(block_variant_t const &block) {
  return std::visit([&](auto &&blk) { return blk.dim(); }, block);
}

int64_t size(block_variant_t const &block) {
  return std::visit([&](auto &&blk) { return blk.size(); }, block);
}

int64_t n_sites(block_variant_t const &block) {
  return std::visit([&](auto &&blk) { return blk.n_sites(); }, block);
}

bool isreal(block_variant_t const &block) {
  return std::visit([&](auto &&blk) { return blk.isreal(); }, block);
}

bool isdistributed(block_variant_t const &block) {
  return std::visit(overload{
                        [&](Spinhalf const &) -> bool { return false; },
                        [&](tJ const &) -> bool { return false; },
                        [&](Electron const &) -> bool { return false; },
#ifdef XDIAG_USE_MPI
                        [&](tJDistributed const &) -> bool { return true; },
#endif
                        [&](auto &&) -> bool { return false; },
                    },
                    block);
}

} // namespace xdiag
