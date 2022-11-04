#include "indexing_variants.h"

namespace hydra::indexing {

idx_t size(SpinhalfIndexing const &idxing) {
  return std::visit(
      variant::overloaded{[&](auto const &idx) -> idx_t { return idx.size(); }},
      idxing);
}

idx_t size(ElectronIndexing const &idxing) {
  return std::visit(
      variant::overloaded{[&](auto const &idx) -> idx_t { return idx.size(); }},
      idxing);
}

idx_t size(tJIndexing const &idxing) {
  return std::visit(
      variant::overloaded{[&](auto const &idx) -> idx_t { return idx.size(); }},
      idxing);
}

} // namespace hydra::indexing
