#include <xdiag/all.hpp>

using namespace xdiag;
using namespace arma;

int main(int argc, char *argv[]) try {
  MPI_Init(&argc, &argv);
  assert(argc == 2);
  int64_t N = atoi(argv[1]);
  int64_t nup = N / 2 - 1;
  int64_t ndn = N / 2 - 1;
  XDIAG_SHOW(N);

  Log.set_verbosity(3);

  tic();
  auto block = tJDistributed(N, nup, ndn);
  toc("creation");

  XDIAG_SHOW(block);
  OpSum ops;
  ops["T"] = 1.0;
  ops["J"] = 1.0;
  for (int s = 0; s < N; ++s) {
    ops += "T" * Op("Hop", {s, (s + 1) % N});
    ops += "J" * Op("SdotS", {s, (s + 1) % N});
  }

  tic();
  double e0 = eigval0(ops, block, 1e-12, 5);
  toc("MVM");
  MPI_Finalize();
} catch (Error e) {
  error_trace(e);
}
