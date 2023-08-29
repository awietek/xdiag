#pragma once

#include <utility>

#include <hydra/blocks/blocks.h>
#include <hydra/operators/bondlist.h>
#include <hydra/states/state.h>

namespace hydra {

double eigval0(BondList const &bondlist, block_variant_t const &block,
               double precision = 1e-12, int64_t max_iterations = 1000,
               bool force_complex = false, int64_t random_seed = 42);

std::pair<double, State>
eig0(BondList const &bondlist, block_variant_t const &block,
     double precision = 1e-12, int64_t max_iterations = 1000,
     bool force_complex = false, int64_t random_seed = 42);

} // namespace hydra
