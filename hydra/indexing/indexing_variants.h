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

namespace hydra::indexing {

// clang-format off
using SpinhalfIndexing = std::variant<
    spinhalf::IndexingSz<uint16_t>,
    spinhalf::IndexingNoSz<uint16_t>,
    spinhalf::IndexingSymmetricSz<uint16_t>,
    spinhalf::IndexingSymmetricNoSz<uint16_t>,
    spinhalf::IndexingSublattice<uint16_t, 1>,
    spinhalf::IndexingSublattice<uint16_t, 2>,
    spinhalf::IndexingSublattice<uint16_t, 3>,
    spinhalf::IndexingSublattice<uint16_t, 4>,
    spinhalf::IndexingSublattice<uint16_t, 5>,
    spinhalf::IndexingSz<uint32_t>,
    spinhalf::IndexingNoSz<uint32_t>,
    spinhalf::IndexingSymmetricSz<uint32_t>,
    spinhalf::IndexingSymmetricNoSz<uint32_t>,
    spinhalf::IndexingSublattice<uint32_t, 1>,
    spinhalf::IndexingSublattice<uint32_t, 2>,
    spinhalf::IndexingSublattice<uint32_t, 3>,
    spinhalf::IndexingSublattice<uint32_t, 4>,
    spinhalf::IndexingSublattice<uint32_t, 5>,
    spinhalf::IndexingSz<uint64_t>,
    spinhalf::IndexingNoSz<uint64_t>,
    spinhalf::IndexingSymmetricSz<uint64_t>,
    spinhalf::IndexingSymmetricNoSz<uint64_t>,
    spinhalf::IndexingSublattice<uint64_t, 1>,
    spinhalf::IndexingSublattice<uint64_t, 2>,
    spinhalf::IndexingSublattice<uint64_t, 3>,
    spinhalf::IndexingSublattice<uint64_t, 4>,
    spinhalf::IndexingSublattice<uint64_t, 5>>;
// clang-format on

idx_t size(SpinhalfIndexing const &idxing);

// clang-format off
using ElectronIndexing =
  std::variant<electron::IndexingNp<uint16_t>,
	       electron::IndexingNoNp<uint16_t>,
	       electron::IndexingSymmetricNp<uint16_t>,
	       electron::IndexingSymmetricNoNp<uint16_t>,
	       electron::IndexingNp<uint32_t>,
	       electron::IndexingNoNp<uint32_t>,
	       electron::IndexingSymmetricNp<uint32_t>,
	       electron::IndexingSymmetricNoNp<uint32_t>,
	       electron::IndexingNp<uint64_t>,
	       electron::IndexingNoNp<uint64_t>,
	       electron::IndexingSymmetricNp<uint64_t>,
	       electron::IndexingSymmetricNoNp<uint64_t>>;
// clang-format on

idx_t size(ElectronIndexing const &idxing);

// clang-format off
using tJIndexing =
  std::variant<tj::IndexingNp<uint16_t>,
	       tj::IndexingSymmetricNp<uint16_t>,
	       tj::IndexingNp<uint32_t>,
	       tj::IndexingSymmetricNp<uint32_t>,
	       tj::IndexingNp<uint64_t>,
	       tj::IndexingSymmetricNp<uint64_t>>;
// clang-format on

idx_t size(tJIndexing const &idxing);

} // namespace hydra::indexing
