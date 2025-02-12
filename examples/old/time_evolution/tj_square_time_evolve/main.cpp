#include <xdiag/all.hpp>

using namespace xdiag;

void measure_density(int nsites, State const &v) {
  for (int i = 0; i < nsites; ++i) {
    auto sz = innerC(Bond("Ntot", i), v);
    printf("%f ", std::real(sz));
  }
  printf("\n");
}

int main() {

  int L = 3;
  double t = 1.0;
  double J = 10.0;

  int nsites = L * L;
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
      ops += Op("Hop", "T", {site, right});
      ops += Op("tJSzSz", "J", {site, right});
      ops += Op("Hop", "T", {site, top});
      ops += Op("tJSzSz", "J", {site, top});
    }
  }
  ops["T"] = t;
  ops["J"] = J;

  auto block = tJ(nsites, nsites / 2, nsites / 2 - 1);
  ops["J"] = J / 2;
  auto [e0, v] = eig0(ops, block);

  measure_density(nsites, v);

  ops["J"] = J;

  
  // Do the time evolution with a step size tau
  double tau = 0.1;
  for (int i = 0; i < 40; ++i) {
    v = time_evolve(ops, v, tau, precision);
    measure_density(nsites, v);
  }

  return EXIT_SUCCESS;
}
