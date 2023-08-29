#include "blocks.h"

namespace hydra {
int64_t size(block_variant_t const &block) {
  return std::visit(overload{[&](auto &&blk) -> idx_t { return blk.size(); }},
                    block);
}

int64_t n_sites(block_variant_t const &block) {
  return std::visit(
      overload{[&](auto &&blk) -> idx_t { return blk.n_sites(); }}, block);
}

bool isreal(block_variant_t const &block) {
  return std::visit(overload{[&](auto &&blk) -> idx_t { return blk.isreal(); }},
                    block);
}

bool iscomplex(block_variant_t const &block) { return !isreal(block); }
} // namespace hydra
