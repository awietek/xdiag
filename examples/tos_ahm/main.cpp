#include <xdiag/all.hpp>

using namespace xdiag;

int main() try {
  say_hello();
  int Lx = 4;
  int Ly = 4;
  int Nmin = 4;
  int Nmax = 2 * Lx * Ly - 4;
  std::vector<double> Us = {0.0, -2.0, -10.0};

  for (double U : Us) {

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
        ops += "mu" * Op("Ntot", s);     // chemical potential
      }
    }
    ops["t"] = 1;
    ops["mu"] = -U / 2;
    ops["U"] = U;

    for (int N = Nmin; N <= Nmax; N += 2) {
      int nup = N / 2;
      int ndn = N / 2;
      auto block = Electron(Lx * Ly, nup, ndn); // create Hilbert space

      auto res = eigs_lanczos(ops, block, 3);

      std::string outfile = fmt::format(
          "data/tos_ahm/U({})_N({})_Lx({})_Ly({}).h5", U, N, Lx, Ly);
      auto f = FileH5(outfile, "w!");
      f["spectrum"] = res.eigenvalues;
    }
  }
  return 0;
} catch (Error e) {
  error_trace(e);
}
