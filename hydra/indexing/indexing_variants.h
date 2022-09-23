#pragma once
#include <variant>

#include <hydra/indexing/spinhalf/indexing_no_sz.h>
#include <hydra/indexing/spinhalf/indexing_sublattice.h>
#include <hydra/indexing/spinhalf/indexing_symmetric_no_sz.h>
#include <hydra/indexing/spinhalf/indexing_symmetric_sz.h>
#include <hydra/indexing/spinhalf/indexing_sz.h>

#include <hydra/indexing/tj/indexing_np.h>
#include <hydra/indexing/tj/indexing_symmetric_np.h>

#include <hydra/indexing/electron/indexing_no_np.h>
#include <hydra/indexing/electron/indexing_np.h>
#include <hydra/indexing/electron/indexing_symmetric_no_np.h>
#include <hydra/indexing/electron/indexing_symmetric_np.h>

// helper type for visitors
namespace hydra {
template <class... Ts> struct overloaded : Ts... {
  using Ts::operator()...;
};
template <class... Ts> overloaded(Ts...) -> overloaded<Ts...>;

} // namespace hydra

namespace hydra::indexing::spinhalf {

template <typename bit_t>
using Indexing =
    std::variant<IndexingSz<bit_t>, IndexingNoSz<bit_t>,
                 IndexingSymmetricSz<bit_t>, IndexingSymmetricNoSz<bit_t>,
                 IndexingSublattice<bit_t, 1>, IndexingSublattice<bit_t, 2>,
                 IndexingSublattice<bit_t, 3>, IndexingSublattice<bit_t, 4>,
                 IndexingSublattice<bit_t, 5>>;

template <typename bit_t> idx_t size(Indexing<bit_t> const &idxing) {
  return std::visit(
      overloaded{[&](auto const &idx) -> idx_t { return idx.size(); }}, idxing);
}

} // namespace hydra::indexing::spinhalf

namespace hydra::indexing::electron {

template <typename bit_t>
using Indexing =
    std::variant<IndexingNp<bit_t>, IndexingNoNp<bit_t>,
                 IndexingSymmetricNp<bit_t>, IndexingSymmetricNoNp<bit_t>>;

template <typename bit_t> idx_t size(Indexing<bit_t> const &idxing) {
  return std::visit(
      overloaded{[&](auto const &idx) -> idx_t { return idx.size(); }}, idxing);
}

} // namespace hydra::indexing::electron

namespace hydra::indexing::tj {

template <typename bit_t>
using Indexing = std::variant<IndexingNp<bit_t>, IndexingSymmetricNp<bit_t>>;

template <typename bit_t> idx_t size(Indexing<bit_t> const &idxing) {
  return std::visit(
      overloaded{[&](auto const &idx) -> idx_t { return idx.size(); }}, idxing);
}

} // namespace hydra::indexing::tj
