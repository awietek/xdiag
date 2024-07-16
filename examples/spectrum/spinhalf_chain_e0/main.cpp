#include <xdiag/all.hpp>

using namespace xdiag;

int main() try {
  int N = 16;
  int nup = N / 2;
  Spinhalf block(N, nup);

  // Define the nearest-neighbor Heisenberg model
  OpSum ops;
  for (int i = 0; i < N; ++i) {
    ops += Op("HB", "J", {i, (i + 1) % N});
  }
  ops["J"] = 1.0;

  set_verbosity(2);                  // set verbosity for monitoring progress
  double e0 = eigval0(ops, block); // compute ground state energy

  Log("Ground state energy: {:.12f}", e0);

} catch (Error e) {
  error_trace(e);
}
