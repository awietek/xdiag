#include <hydra/all.h>

int main() {
  using namespace hydra;

  lila::Log.set_verbosity(1);

  int n_sites = 7;
  BondList bondlist;
  for (int s = 0; s < n_sites; ++s) {
    bondlist << Bond("HOP", "T", {s, (s + 1) % n_sites});
    bondlist << Bond("ISING", "JZ", {s, (s + 1) % n_sites});
    bondlist << Bond("EXCHANGE", "JXY", {s, (s + 1) % n_sites});
  }

  // cyclic group as space group
  std::vector<std::vector<int>> permutations;
  for (int sym = 0; sym < n_sites; ++sym) {
    std::vector<int> permutation;
    for (int site = 0; site < n_sites; ++site) {
      int newsite = (site + sym) % n_sites;
      permutation.push_back(newsite);
    }
    permutations.push_back(permutation);
  }
  auto space_group = PermutationGroup(permutations);

  // int n_sites = 16;
  // std::string lfile = std::string("square.") + std::to_string(n_sites) +
  //                     std::string(".tJ.fsl.pbc.lat");
  // auto irrep = read_represenation(lfile, "Gamma.D4.A1");
  // auto permutations = read_permutations(lfile);
  // auto space_group = PermutationGroup(permutations);
  // auto bondlist = read_bondlist(lfile);

  Couplings cpls;
  cpls["T"] = 0.0;
  cpls["JZ"] = 0.0;
  cpls["JXY"] = 1.0;

  int n_up = n_sites / 2;
  int n_dn = n_sites - n_up;
  {

    for (int k = 0; k < n_sites; ++k) {
      lila::Log("k={}", k);
      // Create irrep
      std::vector<complex> chis;
      for (int l = 0; l < n_sites; ++l)
        chis.push_back({std::cos(2 * M_PI * l * k / n_sites),
                        std::sin(2 * M_PI * l * k / n_sites)});
      auto irrep = Representation(chis);
      auto block = tJSymmetric(n_sites, n_up, n_dn, space_group, irrep);
      auto H = MatrixCplx(bondlist, cpls, block, block);
      // LilaPrint(H);
      LilaPrint(lila::EigenvaluesSym(H));
      lila::Log("--------------------");


    }
    // lila::tic();
    // double e0 = E0Real(bondlist, cpls, block);
    // lila::toc();

    // lila::Log.out("TJ e0: {}", e0);

    auto block2 = Spinhalf(n_sites, n_up);
    auto H2 = MatrixReal(bondlist, cpls, block2, block2);
    // LilaPrint(H2);
    LilaPrint(lila::EigenvaluesSym(H2));

    // lila::tic();
    // auto block3 = Electron(n_sites, n_up, n_dn);
    // lila::toc("build block");
    // auto H3 = MatrixReal(bondlist, cpls, block3, block3);
    // LilaPrint(H3);
    // LilaPrint(lila::EigenvaluesSym(H3));

    // lila::tic();
    // double e02 = E0Real(bondlist, cpls, block2);
    // lila::toc();

    // lila::Log.out("Spinhalf e0: {}", e02);
  }

  return EXIT_SUCCESS;
}
