#pragma once

#include <variant>

#include <hydra/blocks/electron/electron.h>
#include <hydra/blocks/spinhalf/spinhalf.h>
#include <hydra/blocks/tj/tj.h>
#include <hydra/common.h>

namespace hydra {

using block_variant_t = std::variant<Spinhalf, tJ, Electron>;
int64_t size(block_variant_t const &block);
int64_t n_sites(block_variant_t const &block);
bool isreal(block_variant_t const &block);
bool iscomplex(block_variant_t const &block);

} // namespace hydra
