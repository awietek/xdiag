#include <xdiag/all.hpp>

using namespace xdiag;

int main() try {

  int n_sites = 16;
  int nup = n_sites / 2;

  // Define the Hilbert space block
  auto block = Spinhalf(n_sites, nup);

  // Define the nearest-neighbor Heisenberg Hamiltonian
  BondList bonds;
  for (int i = 0; i < n_sites; ++i) {
    bonds << Bond("HB", "J", {i, (i + 1) % n_sites});
  }

  // Set the coupling constant "J" to one
  bonds["J"] = 1.0;

  // Compute and print the ground state energy
  double e0 = eigval0(bonds, block);
  Log("Ground state energy: {:.12f}", e0);

} catch (std::exception const &e) {
  traceback(e);
}
