#include "../../catch.hpp"

#include <hydra/all.h>

TEST_CASE("kitaev_gamma", "[blocks][spinhalf_symmetric]") {
  using namespace hydra;
  using namespace arma;

  std::string lfile = "data/kitaev_gamma/lattice-files/"
                      "honeycomb.8.HeisenbergKitaevGamma.fsl.lat";
  auto bonds = read_bondlist(lfile);
  auto group = PermutationGroup(read_permutations(lfile));

  // double K = -.30901699437494742412;
  // double G = .95105651629515357210;

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


  // double K = -.30901699437494742412;
  // double G = .95105651629515357210; 
  // std::vector<std::pair<std::string, double>> irrep_names_e0 = {
  //     {"Gamma.C2.A", -3.1765766652975568896},
  //     {"Gamma.C2.B", -2.0558245985778302867},
  //     {"M0.C2.A", -2.3210367200223522843},
  //     {"M0.C2.B", -2.8044982600737728973},
  //     {"M1.C2.A", -2.3210367200223518402},
  //     {"M1.C2.B", -2.8044982600737733414},
  //     {"M2.C2.A", -2.3210367200223536166},
  //     {"M2.C2.B", -2.8044982600737742295}};



  double K = -1.0;
  double G = 0.0;
  std::vector<std::pair<std::string, double>> irrep_names_e0 = {
      {"Gamma.C2.A", -1.732050807568876083},
      {"Gamma.C2.B", -0.85463767971846216209},
      {"M0.C2.A", -1.4516059629557762634},
      {"M0.C2.B", -1.4516059629557771515},
      {"M1.C2.A", -1.4516059629557773736},
      {"M1.C2.B", -1.4516059629557787058},
      {"M2.C2.A", -1.4516059629557762634},
      {"M2.C2.B", -1.4516059629557764854}};
  

  bonds["J"] = 0.;
  bonds["K"] = K;
  bonds["G"] = G;
  for (auto [name, e0_reference] : irrep_names_e0) {
    auto irrep = read_represenation(lfile, name);
    auto block = Spinhalf(8, group, irrep);
    // double e0_computed = e0(bonds, block);
    cx_mat H = matrix(bonds, block);

    vec eigs;
    eig_sym(eigs, H);
    double e0_computed = eigs(0);

    HydraPrint(name);
    HydraPrint(e0_reference);
    HydraPrint(e0_computed);
      
  }
}
