#include <hydra/all.h>

void measure_density(int n_sites, hydra::State const &v) {
  using namespace hydra;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  for (int i = 0; i < n_sites; ++i) {
    auto sz = inner(Bond("NUMBER", i), v);
    if (rank == 0) {
      printf("%.6f ", std::real(sz));
    }
  }
  if (rank == 0) {
    printf("\n");
  }
}

int main() {
  using namespace hydra;

  int L = 4;
  int W = 4;
  double t = 1.0;
  double J = 0.4;
  double mu_0 = 2;

  int n_sites = L * L;
  double precision = 1e-12;

  // Create square lattice t-J model
  BondList bonds;
  for (int x = 0; x < L; ++x) {
    for (int y = 0; y < W; ++y) {
      int nx = (x + 1) % L;
      int ny = (y + 1) % W;

      int site = y * L + x;
      int right = y * L + nx;
      int top = ny * L + x;
      bonds << Bond("HOP", "T", {site, right});
      bonds << Bond("TJISING", "J", {site, right});
      bonds << Bond("HOP", "T", {site, top});
      bonds << Bond("TJISING", "J", {site, top});

      if (x < L / 2) {
        bonds << Bond("NUMBER", "MUPLUS", site);
      } else {
        bonds << Bond("NUMBER", "MUNEG", site);
      }
    }
  }
  bonds["T"] = t;
  bonds["J"] = J;
  bonds["MUPLUS"] = mu_0;
  bonds["MUNEG"] = mu_0;

  auto block = tJDistributed(n_sites, n_sites / 2, n_sites / 2 - 1);
  auto [e0, v] = eig0(bonds, block);

  bonds["MUPLUS"] = 0;
  bonds["MUNEG"] = 0;

  measure_density(n_sites, v);

  // Do the time evolution with a step size tau
  double tau = 0.1;
  for (int i = 0; i < 40; ++i) {
    v = time_evolve(bonds, v, tau, precision);
    measure_density(n_sites, v);
  }

  return EXIT_SUCCESS;
}
