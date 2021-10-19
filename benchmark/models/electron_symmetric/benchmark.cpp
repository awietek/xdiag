#include <hydra/all.h>

int main() {
  using namespace hydra;

  lila::Log.set_verbosity(1);

  // int n_sites = 16;
  // std::string lfile = std::string("triangular.") + std::to_string(n_sites) +
  //                     std::string(".t1t2.pbc.lat");
  // auto irrep = read_represenation(lfile, "Gamma.D6.A1");



  int n_sites = 18;
  std::string lfile = std::string("triangular.") + std::to_string(n_sites) +
                      std::string(".t1t2.pbc.lat");
  auto irrep = read_represenation(lfile, "Gamma.C2.A");


  int n_up = n_sites / 2;
  int n_dn = n_sites - n_up;
  auto bondlist = read_bondlist(lfile);
  Couplings cpls;
  cpls["T1"] = 1.0;
  cpls["U"] = 5.0;
  auto permutations = read_permutations(lfile);
  auto space_group = PermutationGroup(permutations);
  lila::tic();
  auto block = ElectronSymmetric(n_sites, n_up, n_dn, space_group, irrep);
  lila::toc("build block");

  lila::tic();
  double e0 = E0Real(bondlist, cpls, block);
  lila::toc();
  lila::Log.out("e0: {}", e0);

  return EXIT_SUCCESS;
}
