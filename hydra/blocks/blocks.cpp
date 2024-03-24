#include "blocks.h"

#include <hydra/parallel/mpi/cdot_distributed.h>

namespace hydra {

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

bool iscomplex(block_variant_t const &block) { return !isreal(block); }

bool isdistributed(block_variant_t const &block) {
  return std::visit(overload{
                        [&](Spinhalf const &) -> bool { return false; },
                        [&](tJ const &) -> bool { return false; },
                        [&](Electron const &) -> bool { return false; },
#ifdef HYDRA_USE_MPI
                        [&](tJDistributed const &) -> bool { return true; },
#endif
                        [&](auto &&) -> bool { return false; },
                    },
                    block);
}

} // namespace hydra
