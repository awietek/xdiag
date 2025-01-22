#include <xdiag/all.hpp>

void measure_density(int nsites, xdiag::State const &v) {
  using namespace xdiag;

  for (int i = 0; i < nsites; ++i) {
    auto sz = inner(Bond("Ntot", i), v);
    printf("%f ", real(sz));
  }
  printf("\n");
}

int main() {
  using namespace xdiag;  
  int nsites = 8;
  double t = 1.0;
  double U = 8.0;

  double precision = 1e-12;

  OpSum ops;
  for (int i = 0; i < nsites; ++i) {
    ops << Op("Hop", "T", {i, (i + 1) % nsites});
  }
  ops["T"] = t;
  ops["U"] = U;

  auto pstate = ProductState();
  for (int i = 0; i < nsites; ++i) {
    if ((i % 2) == 0) {
      pstate << "Dn";
    } else {
      pstate << "Up";
    }
  }
  pstate[3] = "UpDn";

  auto block = Electron(nsites, 4, 5);
  auto v = State(block, pstate);

  measure_density(nsites, v);

  // Do the time evolution with a step size tau
  double tau = 0.1;
  for (int i = 0; i < 40; ++i) {
    v = time_evolve(ops, v, tau, precision);
    measure_density(nsites, v);
  }

  return EXIT_SUCCESS;
}
