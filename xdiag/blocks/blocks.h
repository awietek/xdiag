#pragma once

#include <functional>
#include <variant>

#include <xdiag/extern/armadillo/armadillo>

#include <xdiag/common.h>

#include <xdiag/blocks/electron/electron.h>
#include <xdiag/blocks/spinhalf/spinhalf.h>
#include <xdiag/blocks/tj/tj.h>
#include <xdiag/blocks/tj_distributed/tj_distributed.h>

namespace xdiag {

#ifdef XDIAG_USE_MPI
using block_variant_t = std::variant<Spinhalf, tJ, Electron, tJDistributed>;
#else
using block_variant_t = std::variant<Spinhalf, tJ, Electron>;
#endif

int64_t dim(block_variant_t const &block);
int64_t size(block_variant_t const &block);
int64_t n_sites(block_variant_t const &block);
bool isreal(block_variant_t const &block);
bool iscomplex(block_variant_t const &block);
bool isdistributed(block_variant_t const &block);


} // namespace xdiag
