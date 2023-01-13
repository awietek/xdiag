#include "blocks.h"
#include <hydra/random/hashes.h>

namespace hydra {

Block::Block(Spinhalf const &variant) : variant_(variant) {}
Block::Block(tJ const &variant) : variant_(variant) {}
Block::Block(Electron const &variant) : variant_(variant) {}
Block::Block(block_variant_t const &variant) : variant_(variant) {}

idx_t Block::size() const {
  return std::visit(
      variant::overloaded{[&](auto &&blk) -> idx_t { return blk.size(); }},
      variant_);
}

uint64_t Block::hash() const {
  return std::visit(variant::overloaded{[&](auto &&blk) -> idx_t {
                      return random::hash(blk);
                    }},
                    variant_);
}

block_variant_t &Block::variant() { return variant_; }
block_variant_t const &Block::variant() const { return variant_; }

bool Block::operator==(Block const &rhs) const {
  return (variant_ == rhs.variant_);
}

bool Block::operator!=(Block const &rhs) const { return !operator==(rhs); }

} // namespace hydra
