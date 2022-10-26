#include "blocks.h"
#include <hydra/random/hashes.h>

namespace hydra {
idx_t size(Block const &block) {
  return std::visit(
      overloaded{[&](auto const &blk) -> idx_t { return blk.size(); }}, block);
}

idx_t hash(Block const &block) {
  return std::visit(
      overloaded{[&](auto const &blk) -> idx_t { return random::hash(block); }},
      block);
}

} // namespace hydra
