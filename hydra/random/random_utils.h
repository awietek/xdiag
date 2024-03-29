#pragma once

#include <random>

#include <hydra/common.h>

#include <hydra/extern/armadillo/armadillo>

namespace hydra::random {

double random_uniform_real(std::mt19937 &gen, double a = 0.0, double b = 1.0);

double boxmuller_single(std::mt19937 &gen);

// Comment: std::normal_distribution is not deterministic with discard
// so hydra uses own boxmuller algorithm to create random normal numbers
double random_normal(std::mt19937 &gen, double mean = 0.0,
                     double variance = 1.0);

int random_uniform_real_discard();
int random_normal_discard();

template <typename coeff_t>
void fill_random_normal_vector(arma::Col<coeff_t> &v, int seed);

} // namespace hydra::random
