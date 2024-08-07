#include <xdiag/all.hpp>

int main() {
  using namespace xdiag;
  using namespace arma;

  int n_sites = 8;

  // Create the Hamiltonian
  OpSum ops;
  auto T_mat = mat(n_sites, n_sites, fill::randn);
  auto J_mat = mat(n_sites, n_sites, fill::randn);
  for (int i = 0; i < n_sites; ++i) {
    for (int j = 0; j < n_sites; ++j) {
      if (i != j) {

        // Define hoppings
        std::string T = fmt::format("T{}{}", i, j);
        ops += Op("HOP", T, {i, j});
        ops[T] = T_mat(i, j);

        // Define Heisenberg interactions
        std::string J = fmt::format("J{}{}", i, j);
        ops += Op("HB", J, {i, j});
        ops[J] = J_mat(i, j);
      }
    }
  }

  // Define block with two holes
  int n_up = n_sites / 2 - 1;
  int n_dn = n_sites / 2 - 1;
  auto block = tJ(n_sites, n_up, n_dn);

  // Compute the full Hamiltonian matrix
  auto H = matrixC(ops, block);

  // Compute all eigenvalues of H
  vec eigs;
  eig_sym(eigs, H);
  eigs.print();

  return EXIT_SUCCESS;
}
