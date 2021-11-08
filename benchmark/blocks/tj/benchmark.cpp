#include <hydra/all.h>

int main() {
  using namespace hydra;

  lila::Log.set_verbosity(1);

  int n_sites = 16;
  std::string lfile = std::string("square.") + std::to_string(n_sites) +
                      std::string(".tJ.fsl.pbc.lat");

  int n_up = n_sites / 2 - 2;
  int n_dn = n_sites / 2 - 2;

  // int n_up = 0; // n_sites / 2;
  // int n_dn = 6;


  auto bondlist = read_bondlist(lfile);
  Couplings cpls;
  cpls["T"] = 1.0;
  cpls["JZ"] = 1.0;
  cpls["JXY"] = 1.0;
  {
    lila::tic();
    auto block = tJ(n_sites, n_up, n_dn);
    lila::toc("build block");



    lila::tic();
    double e0 = E0Real(bondlist, cpls, block);
    lila::toc();

    lila::Log.out("TJ e0: {}", e0);
  }

  return EXIT_SUCCESS;
}
