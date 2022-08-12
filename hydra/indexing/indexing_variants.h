#pragma once
#include <variant>

#include <hydra/indexing/spinhalf/indexing_no_sz.h>
#include <hydra/indexing/spinhalf/indexing_sublattice.h>
#include <hydra/indexing/spinhalf/indexing_symmetric_no_sz.h>
#include <hydra/indexing/spinhalf/indexing_symmetric_sz.h>
#include <hydra/indexing/spinhalf/indexing_sz.h>

#include <hydra/indexing/tj/tj_indexing.h>
#include <hydra/indexing/tj/tj_symmetric_indexing.h>

#include <hydra/indexing/electron/electron_indexing.h>
#include <hydra/indexing/electron/electron_indexing_no_np.h>
#include <hydra/indexing/electron/electron_symmetric_indexing.h>
#include <hydra/indexing/electron/electron_symmetric_indexing_no_np.h>

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

template <typename bit_t> idx_t dimension(Indexing<bit_t> const &idxing) {
  return std::visit(
      overloaded{[&](auto const &idx) -> idx_t { return idx.size(); }}, idxing);
}

} // namespace hydra::indexing::spinhalf
