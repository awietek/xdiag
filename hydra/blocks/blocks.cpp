#include "blocks.h"
#include <hydra/random/hashes.h>

namespace hydra {
idx_t size(Block const &block) {
  return std::visit(
      variant::overloaded{[&](auto &&blk) -> idx_t { return blk.size(); }},
      block);
}

uint64_t hash(Block const &block) {
  return std::visit(variant::overloaded{[&](auto &&blk) -> idx_t {
                      return random::hash(blk);
                    }},
                    block);
}

} // namespace hydra
