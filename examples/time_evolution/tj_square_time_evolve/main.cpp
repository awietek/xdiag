#include <xdiag/all.hpp>

using namespace xdiag;

void measure_density(int n_sites, State const &v) {
  for (int i = 0; i < n_sites; ++i) {
    auto sz = innerC(Bond("NUMBER", i), v);
    printf("%f ", std::real(sz));
  }
  printf("\n");
}

int main() {

  int L = 3;
  double t = 1.0;
  double J = 10.0;

  int n_sites = L * L;
  double precision = 1e-12;

  // Create square lattice t-J model
  OpSum ops;
  for (int x = 0; x < L; ++x) {
    for (int y = 0; y < L; ++y) {
      int nx = (x + 1) % L;
      int ny = (y + 1) % L;

      int site = y * L + x;
      int right = y * L + nx;
      int top = ny * L + x;
      ops += Op("HOP", "T", {site, right});
      ops += Op("TJISING", "J", {site, right});
      ops += Op("HOP", "T", {site, top});
      ops += Op("TJISING", "J", {site, top});
    }
  }
  ops["T"] = t;
  ops["J"] = J;

  auto block = tJ(n_sites, n_sites / 2, n_sites / 2 - 1);
  ops["J"] = J / 2;
  auto [e0, v] = eig0(ops, block);

  measure_density(n_sites, v);

  ops["J"] = J;

  
  // Do the time evolution with a step size tau
  double tau = 0.1;
  for (int i = 0; i < 40; ++i) {
    v = time_evolve(ops, v, tau, precision);
    measure_density(n_sites, v);
  }

  return EXIT_SUCCESS;
}
