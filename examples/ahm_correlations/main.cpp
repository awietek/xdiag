#include <xdiag/all.hpp>

using namespace xdiag;

int main() try {
  say_hello();
  double U = -10.0;
  int Lx = 2;
  int Ly = 2;

  // Construct Hamiltonian
  auto ops = OpSum();
  ops += "U" * Op("HubbardU"); // apply Hubbard interaction to all sites

  for (int x = 0; x < Lx; x++) {
    for (int y = 0; y < Ly; y++) {
      int s = x * Ly + y;
      int sx = ((x + 1) % Lx) * Ly + y;
      int sy = x * Ly + (y + 1) % Ly;

      ops += "t" * Op("Hop", {s, sx}); // hopping in x-direction
      ops += "t" * Op("Hop", {s, sy}); // hopping in y-direction
    }
  }
  ops["t"] = 1;
  ops["U"] = U;

  // Create Hilbert space
  int nup = Lx * Ly / 2;
  int ndn = Lx * Ly / 2;
  auto block = Electron(Lx * Ly, nup, ndn);

  // Get ground state
  auto [e0, psi0] = eig0(ops, block);

  // build correlation matrix

  arma::mat corr = arma::mat(Lx, Ly);
  for (int x = 0; x < Lx; x++) {
    for (int y = 0; y < Ly; y++) {
      int s = x * Ly + y;
      if (s == 1) {
        continue;
      }
      auto op = Op("NtotNtot", {1, s});
      corr(s) = inner(op, psi0);
    }
  }

  return 0;
} catch (Error e) {
  error_trace(e);
}
