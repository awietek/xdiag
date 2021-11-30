#include <hydra/all.h>

int main() {
  using namespace hydra;

  lila::Log.set_verbosity(1);

  int n_sites = 20;
  std::string lfile = "square.20.J1J2.fsl.pbc.lat";
  auto permutations = read_permutations(lfile);
  auto space_group = PermutationGroup(permutations);
  auto bondlist = read_bondlist(lfile);
  auto irrep = read_represenation(lfile, "Gamma.C4.A");

  Couplings cpls;
  cpls["T"] = 1.0;
  cpls["J"] = 0.5;

  int nup = n_sites / 2 - 1;
  int ndn = n_sites / 2 -1;

  auto t1 = lila::rightnow();
  auto block = tJSymmetric(n_sites, nup, ndn, space_group, irrep);
  lila::timing(t1, lila::rightnow(), "build block");
  lila::Log("size: {}", block.size());

  lila::tic();
  double e0 = E0Real(bondlist, cpls, block);
  lila::toc();
  lila::Log.out("e0: {}", e0);

  return EXIT_SUCCESS;
}
