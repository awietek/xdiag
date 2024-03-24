#pragma once

#include <functional>
#include <variant>

#include <extern/armadillo/armadillo>

#include <hydra/common.h>

#include <hydra/blocks/electron/electron.h>
#include <hydra/blocks/spinhalf/spinhalf.h>
#include <hydra/blocks/tj/tj.h>
#include <hydra/blocks/tj_distributed/tj_distributed.h>

namespace hydra {

#ifdef HYDRA_USE_MPI
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


} // namespace hydra
