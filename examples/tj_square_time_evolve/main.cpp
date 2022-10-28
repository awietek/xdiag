#include <hydra/all.h>

void measure_density(int n_sites, hydra::StateCplx const &v) {
  using namespace hydra;

  for (int i = 0; i < n_sites; ++i) {
    auto sz = inner(Bond("NUMBER", i), v);
    printf("%f ", std::real(sz));
  }
  printf("\n");
}

int main() {
  using namespace hydra;

  int L = 5;
  double t = 1.0;
  double J = 0.4;

  int n_sites = L * L;
  double precision = 1e-12;

  // Create square lattice t-J model
  BondList bonds;
  for (int x = 0; x < L; ++x) {
    for (int y = 0; y < L; ++y) {
      int nx = (x + 1) % L;
      int ny = (y + 1) % L;

      int site = y * L + x;
      int right = y * L + nx;
      int top = ny * L + x;
      bonds << Bond("HOP", "T", {site, right});
      bonds << Bond("TJHB", "J", {site, right});
      bonds << Bond("HOP", "T", {site, top});
      bonds << Bond("TJHB", "J", {site, top});
    }
  }
  bonds["T"] = t;
  bonds["J"] = J;

  // Create initial state
  auto pstate = ProductState();
  for (int i = 0; i < n_sites; ++i) {
    if ((i % 2) == 0) {
      pstate << "Up";
    } else {
      pstate << "Dn";
    }
  }
  pstate[12] = "Emp";
  HydraPrint(pstate);

  auto block = tJ(n_sites, 12, 12);
  auto v = State(block, pstate);
  
  measure_density(n_sites, v);

  
  // Do the time evolution with a step size tau
  double tau = 0.1;
  for (int i = 0; i < 100; ++i) {
    v = time_evolve(bonds, v, tau, precision);
    measure_density(n_sites, v);
  }

  return EXIT_SUCCESS;
}
