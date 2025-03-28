#include <xdiag/all.hpp>

using namespace xdiag;
using namespace arma;

int main(int argc, char *argv[]) try {
  MPI_Init(&argc, &argv);
  int N = atoi(argv[1]);
  int nup = N / 2;
  XDIAG_SHOW(N);

  Log.set_verbosity(3);

  tic();
  auto block = SpinhalfDistributed(N, nup);
  toc("creation");

  XDIAG_SHOW(block);
  OpSum ops;
  for (int i = 0; i < N; ++i) {
    ops += Op("SdotS", {i, (i + 1) % N});
  }

  tic();
  double e0 = eigval0(ops, block, 1e-12, 5);
  toc("MVM");
  MPI_Finalize();
} catch (Error e) {
  error_trace(e);
}
