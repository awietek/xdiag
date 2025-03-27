#include <xdiag/all.hpp>

using namespace xdiag;
using namespace arma;

int main(int argc, char *argv[]) try {
  int N = 13;
  auto block = Spinhalf(N);
  XDIAG_SHOW(block);
  OpSum ops;
  for (int i = 0; i < N; ++i) {
    ops += Op("SdotS", {i, (i + 1) % N});
  }
  tic();
  auto H = matrix(ops, block);
  toc("creation");

  vec eigval;
  mat eigvec;
  tic();
  eig_sym(eigval, eigvec, H);
  toc("diagonalization");

} catch (Error e) {
  error_trace(e);
}
