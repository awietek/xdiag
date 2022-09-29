#include "../../catch.hpp"

#include <hydra/all.h>

TEST_CASE("kitaev_gamma", "[blocks][spinhalf_symmetric]") {
  using namespace hydra;
  using namespace arma;

  std::string lfile = "data/kitaev_gamma/lattice-files/"
                      "honeycomb.8.HeisenbergKitaevGamma.fsl.lat";
  auto bonds = read_bondlist(lfile);
  auto group = PermutationGroup(read_permutations(lfile));

  double K = .30901699437494742412;
  double G = .95105651629515357210;

  bonds["J"] = 0.;
  bonds["K"] = K;
  bonds["G"] = G;

  cx_mat sx(mat({{0., 0.5}, {0.5, 0.}}), mat({{0., 0.}, {0., 0.}}));
  cx_mat sy(mat({{0., 0.}, {0., 0.}}), mat({{0., -0.5}, {0.5, 0.}}));
  cx_mat sz(mat({{0.5, 0.0}, {0.0, 0.5}}), mat({{0., 0.}, {0., 0.0}}));

  cx_mat kx = kron(sx, sx);
  cx_mat ky = kron(sy, sy);
  cx_mat kz = kron(sz, sz);
  cx_mat gx = kron(sy, sz) + kron(sz, sy);
  cx_mat gy = kron(sx, sz) + kron(sz, sx);
  cx_mat gz = kron(sx, sy) + kron(sy, sx);
  
  bonds.set_matrix("GENERICKITAEVX", kx);
  bonds.set_matrix("GENERICKITAEVY", ky);
  bonds.set_matrix("GENERICKITAEVZ", kz);
  bonds.set_matrix("GENERICGAMMAX", gx);
  bonds.set_matrix("GENERICGAMMAY", gy);
  bonds.set_matrix("GENERICGAMMAZ", gz);
  
  std::vector<std::pair<std::string, double>> irrep_names_e0 = {
      {"Gamma.C2.A", -3.1765766652975568896},
      {"Gamma.C2.B", -2.0558245985778302867},
      {"M0.C2.A", -2.3210367200223522843},
      {"M0.C2.B", -2.8044982600737728973},
      {"M1.C2.A", -2.3210367200223518402},
      {"M1.C2.B", -2.8044982600737733414},
      {"M2.C2.A", -2.3210367200223536166},
      {"M2.C2.B", -2.8044982600737742295}};

  for (auto [name, e0_reference] : irrep_names_e0) {
    auto irrep = read_represenation(lfile, name);
    auto block = Spinhalf(8, group, irrep);
    double e0_computed = e0(bonds, block);
    HydraPrint(block);
    HydraPrint(e0_reference);
    HydraPrint(e0_computed);
      
  }
}
