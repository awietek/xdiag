#include <hydra/all.h>

void measure_density(int n_sites, hydra::StateCplx const &v) {
  using namespace hydra;

  for (int i = 0; i < n_sites; ++i) {
    auto sz = inner(Bond("NUMBER", i), v);
    printf("%f ", real(sz));
  }
  printf("\n");
}

int main() {
  using namespace hydra;

  int n_sites = 10;
  double t = 1.0;
  double U = 8.0;

  double precision = 1e-12;

  BondList bonds;
  for (int i = 0; i < n_sites; ++i) {
    bonds << Bond("HOP", "T", {i, (i + 1) % n_sites});
  }
  bonds["T"] = t;
  bonds["U"] = U;

  auto pstate = ProductState();
  for (int i = 0; i < n_sites; ++i) {
    if ((i % 2) == 0) {
      pstate << "Dn";
    } else {
      pstate << "Up";
    }
  }
  pstate[5] = "UpDn";

  auto block = Electron(n_sites, 5, 6);
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
