#include <hydra/all.h>

int main() {
  using namespace hydra;

  lila::Log.set_verbosity(1);


  int n_sites = 16;
  std::string lfile = std::string("square.") + std::to_string(n_sites) +
                      std::string(".tJ.fsl.pbc.lat");
  auto irrep = read_represenation(lfile, "Gamma.D4.A1");


  int n_up = n_sites / 2 - 1;
  int n_dn = n_sites / 2 - 1;
  auto bondlist = read_bondlist(lfile);
  Couplings cpls;
  cpls["T"] = 1.0;
  cpls["J"] = 1.0;
  auto permutations = read_permutations(lfile);
  auto space_group = PermutationGroup(permutations);
  lila::tic();
  auto block = tJ(n_sites, n_up, n_dn);
  lila::toc("build block");

  lila::tic();
  double e0 = E0Real(bondlist, cpls, block);
  lila::toc();
  lila::Log.out("e0: {}", e0);

  return EXIT_SUCCESS;
}
