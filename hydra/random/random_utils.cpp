#include "random_utils.h"
#include <math.h>

#ifdef _OPENMP
#include <hydra/parallel/omp/omp_utils.h>
#include <hydra/random/hash_functions.h>
#endif

namespace hydra::random {

double random_uniform_real(std::mt19937 &gen, double a, double b) {
  return std::uniform_real_distribution<double>(a, b)(gen);
}

double boxmuller_single(std::mt19937 &gen) {
  constexpr double twopi = 2.0 * 3.14159265358979323846;

  double u1 = 1.0 - random_uniform_real(gen); // [0, 1) -> (0, 1]
  double u2 = random_uniform_real(gen);
  double radius = std::sqrt(-2 * std::log(u1));
  double theta = twopi * u2;
  return radius * std::cos(theta);
  // return radius * std::sin(theta);
}

// Comment: std::normal_distribution is not deterministic with discard
// so hydra uses own boxmuller algorithm to create random normal numbers
double random_normal(std::mt19937 &gen, double mean, double variance) {
  double z = boxmuller_single(gen);
  return mean + variance * z;
}

template <class distro_f> int random_discard(distro_f distro) {
  std::mt19937 gen(42);
  double r1 = distro(gen);
  double r2 = distro(gen);
  (void)r1;
  // Log("r1: {} r2: {}", r1, r2);

  for (int i = 0; i < 10; ++i) {
    std::mt19937 gen2(42);
    gen2.discard(i);
    double r = distro(gen2);
    // Log("r: {}", r);

    if (std::abs(r - r2) < 1e-12) {
      return i;
    }
  }
  Log.err("Error initializing RNG: could not determine "
          "random_discard");
  return 0;
}

int random_uniform_real_discard() {
  return random_discard(
      [](std::mt19937 &gen) { return random_uniform_real(gen); });
}

int random_normal_discard() {
  return random_discard([](std::mt19937 &gen) { return random_normal(gen); });
}

template <typename coeff_t>
void fill_random_normal_vector(arma::Col<coeff_t> &v, int seed) {
#ifdef _OPENMP
#pragma omp parallel
  {

    auto [start, end] = omp::get_omp_start_end(v.size());
    int discard = random::random_normal_discard();

    if constexpr (isreal<coeff_t>()) {
      std::mt19937 gen(seed);
      gen.discard(discard * start);

      for (idx_t idx = start; idx < end; ++idx) {
        v(idx) = random::random_normal(gen);
      }

    } else {
      std::mt19937 genr(seed);
      std::mt19937 geni(random::hash_fnv1((uint32_t)seed));
      genr.discard(discard * start);
      geni.discard(discard * start);

      for (idx_t idx = start; idx < end; ++idx) {
        auto r = random::random_normal(genr);
        auto i = random::random_normal(geni);
        v(idx) = {r, i};
      }
    }
  }
#else
  std::mt19937 gen(seed);
  for (arma::uword idx = 0; idx < v.size(); ++idx) {
    v(idx) = random::random_normal(gen);
  }
#endif
}

template void fill_random_normal_vector<double>(arma::Col<double> &v, int seed);
template void fill_random_normal_vector<complex>(arma::Col<complex> &v,
                                                 int seed);

} // namespace hydra::random
