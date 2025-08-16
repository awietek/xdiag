#include <xdiag/all.hpp>

using namespace xdiag;

int main(int argc, char *argv[]) try {
  assert(argc == 2);
  int64_t nsites = atoi(argv[1]);
  int64_t nup = nsites / 2;

  // say_hello();
  set_verbosity(2);

  // build Hamiltonian
  auto ops = OpSum();
  for (int i = 0; i < nsites; i++) {
    ops += "J" * Op("SdotS", {i, (i + 1) % nsites});
  }
  ops["J"] = 1.0;

  // get group
  std::vector<int64_t> T_perm(nsites, 0);
  for (int i = 0; i < nsites - 1; i++) {
    T_perm[i] = i + 1;
  };
  auto T = Permutation(T_perm);

  // step 2: define cyclic group
  std::vector<Permutation> perms(nsites);
  for (int i = 0; i < nsites; i++) {
    perms[i] = pow(T, i);
  };
  auto group = PermutationGroup(perms);
  auto irrep = Representation(group);

  tic();
  auto block = Spinhalf(nsites, nup, irrep);
  // auto block = Spinhalf(nsites);
  toc("Block creation");

  XDIAG_SHOW(block);

  // tic();
  // auto coo = coo_matrix(ops, block);
  // toc("COO matrix creation");

  tic();
  auto csr = csr_matrix(ops, block);
  toc("CSR matrix creation");

  // tic();
  // auto csc = csc_matrix(ops, block);
  // toc("CSC matrix creation");

  auto x = arma::vec(csr.ncols, arma::fill::randu);
  auto y = arma::vec(csr.nrows, arma::fill::zeros);
  tic();
  apply(csr, x, y);
  toc("CSR matrix multiply");
  
  tic();
  double e0 = eigval0(ops, block, 1e-12, 1);
  toc("MVM");

} catch (Error e) {
  error_trace(e);
}
