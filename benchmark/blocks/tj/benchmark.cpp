#include <hydra/all.h>

int main() {
  using namespace hydra;

  lila::Log.set_verbosity(1);

  int n_sites = 9;
  std::string lfile = std::string("square.") + std::to_string(n_sites) +
                      std::string(".tJ.fsl.pbc.lat");

  int n_up = 6; // n_sites / 2;
  // int n_dn = 5; // n_sites - n_up;
  int n_dn = n_sites - n_up;

  auto bondlist = read_bondlist(lfile);
  Couplings cpls;
  cpls["T"] = 0.0;
  cpls["JZ"] = 0.0;
  cpls["JXY"] = 1.0;
  {
    lila::tic();
    auto block = tJ(n_sites, n_up, n_dn);
    lila::toc("build block");

    lila::tic();
    auto H = MatrixReal(bondlist, cpls, block, block);
    lila::toc();

    lila::tic();
    auto eigs = lila::EigenvaluesSym(H);
    lila::toc();

    lila::Log.out("TJ e0: {}", eigs(0));
  }

  {
    lila::tic();
    auto block = Spinhalf(n_sites, n_up);
    lila::toc("build block");

    lila::tic();
    auto H = MatrixReal(bondlist, cpls, block, block);
    lila::toc();

    lila::tic();
    auto eigs = lila::EigenvaluesSym(H);
    lila::toc();

    lila::Log.out("Spinhalf e0: {}", eigs(0));
  }

  // {
  //   lila::tic();
  //   auto block = Electron(n_sites, n_up, n_dn);
  //   lila::toc("build block");

  //   lila::tic();
  //   auto H = MatrixReal(bondlist, cpls, block, block);
  //   lila::toc();

  //   lila::tic();
  //   auto eigs = lila::EigenvaluesSym(H);
  //   lila::toc();

  //   lila::Log.out("Electron e0: {}", eigs(0));
  // }

  return EXIT_SUCCESS;
}
