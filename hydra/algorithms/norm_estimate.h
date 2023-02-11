#pragma once

#include <functional>

#include <hydra/blocks/blocks.h>
#include <hydra/common.h>
#include <hydra/operators/bondlist.h>

namespace hydra {

double norm_estimate(std::function<arma::vec(arma::vec const &)> const &apply_A,
                     idx_t N, int n_max_attempts = 100);

} // namespace hydra
