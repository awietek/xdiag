#include <xdiag/all.hpp>

using namespace xdiag;

void measure_density(int n_sites, State const &v) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  for (int i = 0; i < n_sites; ++i) {
    complex sz = innerC(Bond("NUMBER", i), v);
    if (rank == 0) {
      printf("%.6f ", std::real(sz));
    }
  }
  if (rank == 0) {
    printf("\n");
  }
}

int main(int argc, char **argv) try {
  MPI_Init(&argc, &argv);

  int L = 6;
  int W = 4;
  double t = 1.0;
  double J = 0.1;
  double mu_0 = 10;

  int n_sites = L * W;
  double precision = 1e-12;

  // Create square lattice t-J model
  BondList bonds;
  for (int x = 0; x < L-1; ++x) {
    for (int y = 0; y < W; ++y) {
      int nx = (x + 1) % L;
      int ny = (y + 1) % W;

      int site = x * W + y;
      int right = nx * W + y;
      int top = x * W + ny;
      bonds << Bond("HOP", "T", {site, right});
      bonds << Bond("EXCHANGE", "J", {site, right});
      bonds << Bond("HOP", "T", {site, top});
      bonds << Bond("EXCHANGE", "J", {site, top});


      
      if (x < L / 2) {
	Log("x {} y {} site {} t {} r {} +", x, y, site, top, right);
        bonds << Bond("NUMBER", "MUPLUS", site);
      } else {
	Log("x {} y {} site {} t {} r {} -", x, y, site, top, right);
        bonds << Bond("NUMBER", "MUNEG", site);
      }
    }
  }
  bonds["T"] = t;
  bonds["J"] = J;
  bonds["MUPLUS"] = mu_0;
  bonds["MUNEG"] = mu_0;

  auto block = tJDistributed(n_sites, n_sites / 2 - 1, n_sites / 2 - 1);

  XDIAG_PRINT(block);

  Log.set_verbosity(2);
  tic();
  auto [e0, v] = eig0(bonds, block);
  toc("gs");

  bonds["MUPLUS"] = 0;
  bonds["MUNEG"] = 0;

  measure_density(n_sites, v);

  // Do the time evolution with a step size tau
  double tau = 0.1;
  for (int i = 0; i < 40; ++i) {
    tic();
    v = time_evolve(bonds, v, tau, precision);
    toc("time evolve");
    tic();
    measure_density(n_sites, v);
    toc("measure");
  }

  MPI_Finalize();
  return EXIT_SUCCESS;
} catch (std::exception const &e) {
  traceback(e);
}
