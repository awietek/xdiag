#pragma once
#include <variant>

#include <hydra/indexing/spinhalf/spinhalf_indexing_no_sz.h>
#include <hydra/indexing/spinhalf/spinhalf_indexing_sz.h>
#include <hydra/indexing/spinhalf/spinhalf_symmetric_indexing_no_sz.h>
#include <hydra/indexing/spinhalf/spinhalf_symmetric_indexing_sz.h>

#include <hydra/indexing/tj/tj_indexing.h>
#include <hydra/indexing/tj/tj_symmetric_indexing.h>

#include <hydra/indexing/electron/electron_indexing.h>
#include <hydra/indexing/electron/electron_indexing_no_np.h>
#include <hydra/indexing/electron/electron_symmetric_indexing.h>
#include <hydra/indexing/electron/electron_symmetric_indexing_no_np.h>

namespace hydra::indexing {

template <typename bit_t>
using SpinhalfIndexing =
    std::variant<SpinhalfIndexingSz<bit_t>, SpinhalfIndexingNoSz<bit_t>,
                 SpinhalfSymmetricIndexingSz<bit_t>,
                 SpinhalfSymmetricIndexingNoSz<bit_t>>;

}
