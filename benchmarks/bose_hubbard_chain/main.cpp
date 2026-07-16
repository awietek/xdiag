#include <xdiag/all.hpp>

using namespace xdiag;

int main(int argc, char *argv[]) try {
  assert(argc == 2);
  int64_t nsites = atoi(argv[1]);
  int64_t d = 4;
  int64_t np = nsites * d / 2;

  say_hello();
  set_verbosity(2);

  // Create Hamiltonian
  OpSum ops;
  ops["T"] = 1.0;
  ops["U"] = 1.0;
  for (int s = 0; s < nsites; ++s) {
    ops += "T" * Op("Hop", {s, (s + 1) % nsites});
  }
  ops += "U" * Op("HubbardU");
  auto irrep = cyclic_group_irrep(nsites, 0);

  tic();
  auto block = Boson(nsites, d, np, irrep);
  toc("Block creation");

  XDIAG_SHOW(block);

  tic();
  double e0 = eigval0(ops, block, 1e-12, 20);
  toc("MVM");

  tic();
  auto es = eigvals(ops, block, 3,  1e-12, 20);
  toc("LOBPCG");
  
} catch (Error e) {
  error_trace(e);
}
