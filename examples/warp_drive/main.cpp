#include <hydra/all.h>

void measure_szs(int n_sites, hydra::StateCplx const &v) {
  using namespace hydra;

  for (int i = 0; i < n_sites; ++i) {
    auto sz = inner(Bond("SZ", i), v);
    printf("%f ", std::real(sz));
  }
  printf("\n");
}

int main() {
  using namespace hydra;
  using namespace arma;

  int N = 17;
  double J = 1.0;
  double JW = 0.0;
  int warp_distance = 3;
  double H = 1.0;

  double precision = 1e-12;

  BondList bonds;
  for (int i = 0; i < N - 1; ++i) {
    bonds << Bond("HB", J, {i, (i + 1) % N});
  }
  
  BondList bonds_warped;
  for (int i = 0; i < N - 1; ++i) {
    if ((i % warp_distance) == 0) {
      // Log.warn("warp relay station at: {}", i);
      bonds_warped << Bond("HB", JW, {i, (i + 1) % N});
    } else {
      bonds_warped << Bond("HB", J, {i, (i + 1) % N});
    }
  }

  // // Create Hilbertspace without particle number conservation
  // int nup = N / 2;
  // auto block_gs = Spinhalf(N, nup);
  // auto gs = groundstate(bonds, block_gs);
  // measure_szs(N, gs);

  

  // auto block_exc = Spinhalf(N, nup + 1);
  // auto v = State(block_exc);
  // apply(Bond("S+", 8), gs, v);

  auto block = Spinhalf(N, 1);
  auto pstate = ProductState();
  for (int i = 0; i < N; ++i) {
    pstate << "Dn";
  }
  pstate[8] = "Up";
  // HydraPrint(pstate);
  auto v = StateCplx(block, pstate);
  
  
  measure_szs(N, v);

  // Do the time evolution with a step size tau
  double tau = 0.1;
  for (int i = 0; i < 100; ++i) {
    v = time_evolve(bonds_warped, v, tau, precision);
    measure_szs(N, v);
  }

  return EXIT_SUCCESS;
}
