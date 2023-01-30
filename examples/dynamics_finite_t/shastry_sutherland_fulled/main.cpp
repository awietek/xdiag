#include <hydra/all.h>

int main(int argc, char **argv) {
  using namespace hydra;
  auto lfile = FileToml("shastry.16.HB.J.Jd.fsl.toml", 'r');

  std::string kname = "X0.C1.A";

  auto coords = lfile["Coordinates"].as<arma::mat>();
  auto q = lfile[kname + std::string(".momentum")].as<arma::vec>();

  auto bonds = BondList(lfile["Interactions"]);
  auto group = PermutationGroup(lfile["Symmetries"]);
  auto irrep = Representation(lfile[kname]);

  for (arma::uword i = 0; i < coords.n_rows; ++i) {
    auto r = coords.row(i);
    r.print();
  }

  HydraPrint(bonds);
  HydraPrint(group);
  HydraPrint(irrep);
  HydraPrint(q);

  return EXIT_SUCCESS;
}
