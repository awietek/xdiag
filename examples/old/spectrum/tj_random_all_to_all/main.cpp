#include <xdiag/all.hpp>

int main() {
  using namespace xdiag;
  using namespace arma;

  int nsites = 8;

  // Create the Hamiltonian
  OpSum ops;
  auto T_mat = mat(nsites, nsites, fill::randn);
  auto J_mat = mat(nsites, nsites, fill::randn);
  for (int i = 0; i < nsites; ++i) {
    for (int j = 0; j < nsites; ++j) {
      if (i != j) {

        // Define hoppings
        std::string T = fmt::format("T{}{}", i, j);
        ops += Op("Hop", T, {i, j});
        ops[T] = T_mat(i, j);

        // Define Heisenberg interactions
        std::string J = fmt::format("J{}{}", i, j);
        ops += Op("SdotS", J, {i, j});
        ops[J] = J_mat(i, j);
      }
    }
  }

  // Define block with two holes
  int nup = nsites / 2 - 1;
  int ndn = nsites / 2 - 1;
  auto block = tJ(nsites, nup, ndn);

  // Compute the full Hamiltonian matrix
  auto H = matrixC(ops, block);

  // Compute all eigenvalues of H
  vec eigs;
  eig_sym(eigs, H);
  eigs.print();

  return EXIT_SUCCESS;
}
