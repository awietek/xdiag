#include <xdiag/all.hpp>

int main() {
  using namespace xdiag;
  using namespace arma;

  // Define the number of lattice sites, number of upspins
  int nsites = 8;
  int nup = 4;

  // Create the corresponding block of the Spinhalf Hilbertspace
  auto block = Spinhalf(nsites, nup);

  // Define Heisenberg chain with interaction-type HB and coupling J
  OpSum ops;
  for (int s = 0; s < nsites; ++s) {
    ops << Op("HB", "J", {s, (s + 1) % nsites});
  }

  // Set the coupling J to 1.0
  ops["J"] = 1.0;

  // Compute the full Hamiltonian matrix (real, double precision)
  auto H = matrix(ops, block);

  // Perform a full diagonalization
  vec eigval;
  eig_sym(eigval, H);

  eigval.print();

  return EXIT_SUCCESS;
}
