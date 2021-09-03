#include <hydra/all.h>

int main() {
  using namespace hydra;

  lila::Log.set_verbosity(1);

  // int n_sites = 20;
  // std::string lfile = std::string("square.") + std::to_string(n_sites) +
  //                     std::string(".J1J2.fsl.pbc.lat");
  // auto irrep = read_represenation(lfile, "Gamma.C4.A");

  int n_sites = 32;
  std::string lfile = std::string("square.") + std::to_string(n_sites) +
                      std::string(".J1J2.fsl.pbc.lat");
  auto irrep = read_represenation(lfile, "Gamma.D4.A1");


  int n_up = n_sites / 2;
  auto bondlist = read_bondlist(lfile);
  Couplings cpls;
  cpls["J1"] = 1.0;
  auto permutations = read_permutations(lfile);
  auto space_group = PermutationGroup(permutations);
  lila::tic();
  auto block = SpinhalfSymmetric<uint32, PermutationGroupLookup<uint32>>(n_sites, n_up, space_group, irrep);
  // auto block = SpinhalfSymmetric<uint32, PermutationGroupAction>(n_sites, n_up, space_group, irrep);

  // auto block = SpinhalfSymmetric<uint64, PermutationGroupLookup<uint64>>(n_sites, n_up, space_group, irrep);
  // auto block = SpinhalfSymmetric<uint64, PermutationGroupAction>(n_sites, n_up, space_group, irrep);

  lila::toc("build block");

  lila::tic();
  double e0 = E0Real(bondlist, cpls, block);
  lila::toc();
  lila::Log.out("e0: {}", e0);

  return EXIT_SUCCESS;
}
