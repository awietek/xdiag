#include <xdiag/all.hpp>

using namespace xdiag;

int main(int argc, char *argv[]) try {
  assert(argc == 2);
  int64_t nsites = atoi(argv[1]);
  int64_t nup = nsites / 2;
  int64_t ndn = nsites / 2;

  say_hello();
  set_verbosity(2);

  // Create Hamiltonian
  OpSum ops;
  ops["T"] = 1.0;
  for (int s = 0; s < nsites; ++s) {
    ops += "T" * Op("Hop", {s, (s + 1) % nsites});
  }
  ops += 1.0 * Op("HubbardU");

  // Create syms
  std::vector<Permutation> permutation_array;
  for (int64_t sym = 0; sym < nsites; ++sym) {

    std::vector<int64_t> pv;
    for (int64_t site = 0; site < nsites; ++site) {
      int64_t newsite = (site + sym) % nsites;
      pv.push_back(newsite);
    }
    permutation_array.push_back(Permutation(pv));
  }
  auto group = PermutationGroup(permutation_array);
  auto irrep = Representation(group);

  tic();
  // auto block = Spinhalf(nsites, nup, irrep);
  auto block = Electron(nsites, nup, ndn, irrep);
  toc("Block creation");

  XDIAG_SHOW(block);

  tic();
  double e0 = eigval0(ops, block, 1e-12, 5);
  toc("MVM");

} catch (Error e) {
  error_trace(e);
}
