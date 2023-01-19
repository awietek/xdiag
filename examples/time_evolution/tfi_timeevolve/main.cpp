#include <hydra/all.h>

int main() {
  using namespace hydra;
  using namespace arma;

  int N = 16;
  double J = 1.0;
  double H = 1.0;

  double precision = 1e-12;

  // Create the transverse field ising Hamiltonian
  BondList bonds;
  for (int i = 0; i < N - 1; ++i) {
    bonds << Bond("ISING", J, {i, (i + 1) % N});
  }
  cx_mat sx(mat({ {0., 0.5}, {0.5, 0.} }), mat({ {0., 0.}, {0., 0.} }));
  for (int i = 0; i < N; ++i) {
    bonds << Bond(sx, H, i);
  }

  // Create Hilbertspace without particle number conservation
  auto block = Spinhalf(N);

  // Create all-up staring state
  auto v = zero_state(block);
  v.vector()(block.size() - 1) = 1.0;

  // Create magnetization operator
  BondList mag;
  for (int i = 0; i < N; ++i) {
    mag << Bond("SZ", i);
  }

  // Do the time evolution with a step size tau
  double tau = 0.1;
  for (int i = 0; i < 100; ++i) {
    v = time_evolve(bonds, v, tau, precision);
    complex m = inner(mag, v);
    HydraPrint(m);
  }

  return EXIT_SUCCESS;
}
