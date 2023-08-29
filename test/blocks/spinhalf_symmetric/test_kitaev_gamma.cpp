#include "../../catch.hpp"

#include <hydra/blocks/spinhalf/spinhalf_apply.h>
#include <hydra/algebra/algebra.h>
#include <hydra/algebra/matrix.h>
#include <hydra/algorithms/sparse_diag.h>
#include <hydra/utils/close.h>

void run_kitaev_gamma_test(
    double K, double G,
    std::vector<std::pair<std::string, double>> irrep_names_e0) {
  using namespace hydra;
  using namespace arma;

  std::string lfile =
      HYDRA_DIRECTORY "/misc/data/kitaev_gamma/lattice-files/"
                      "honeycomb.8.HeisenbergKitaevGamma.fsl.lat";
  auto bonds = read_bondlist(lfile);
  auto group = PermutationGroup(read_permutations(lfile));

  cx_mat sx(mat({{0., 0.5}, {0.5, 0.}}), mat({{0., 0.}, {0., 0.}}));
  cx_mat sy(mat({{0., 0.}, {0., 0.}}), mat({{0., -0.5}, {0.5, 0.}}));
  cx_mat sz(mat({{0.5, 0.0}, {0.0, -0.5}}), mat({{0., 0.}, {0., 0.0}}));

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

  bonds["J"] = 0.;
  bonds["K"] = K;
  bonds["G"] = G;
  for (auto [name, e0_reference] : irrep_names_e0) {
    auto irrep = read_representation(lfile, name);
    auto block = Spinhalf(8, group, irrep);
    // double e0_computed = e0(bonds, block);
    cx_mat H = matrixC(bonds, block);

    vec eigs;
    eig_sym(eigs, H);
    double e0_computed = eigs(0);

    // HydraPrint(K);
    // HydraPrint(G);
    // HydraPrint(name);
    // HydraPrint(e0_reference);
    // HydraPrint(e0_computed);
    REQUIRE(close(e0_reference, e0_computed));
  }
}

TEST_CASE("kitaev_gamma", "[blocks][spinhalf_symmetric]") {
  using namespace hydra;

  Log("Testing Kitaev-Gamma model");
  {
    double K = -.30901699437494742412;
    double G = .95105651629515357210;
    std::vector<std::pair<std::string, double>> irrep_names_e0 = {
        {"Gamma.C2.A", -3.1765766652975568896},
        {"Gamma.C2.B", -2.0558245985778302867},
        {"M0.C2.A", -2.3210367200223522843},
        {"M0.C2.B", -2.8044982600737728973},
        {"M1.C2.A", -2.3210367200223518402},
        {"M1.C2.B", -2.8044982600737733414},
        {"M2.C2.A", -2.3210367200223536166},
        {"M2.C2.B", -2.8044982600737742295}};
    run_kitaev_gamma_test(K, G, irrep_names_e0);
  }

  {
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
    run_kitaev_gamma_test(K, G, irrep_names_e0);
  }

  {
    double K = 0.0;
    double G = 1.0;
    std::vector<std::pair<std::string, double>> irrep_names_e0 = {
        {"Gamma.C2.A", -2.9261296954884912225},
        {"Gamma.C2.B", -2.2611516188994311705},
        {"M0.C2.A", -2.3314966606292610862},
        {"M0.C2.B", -2.6582065612437091318},
        {"M1.C2.A", -2.3314966606292633067},
        {"M1.C2.B", -2.6582065612437086877},
        {"M2.C2.A", -2.3314966606292646389},
        {"M2.C2.B", -2.6582065612437109081}};
    run_kitaev_gamma_test(K, G, irrep_names_e0);
  }

  {
    double K = -0.64944804833018365819;
    double G = 0.76040596560003093085;
    std::vector<std::pair<std::string, double>> irrep_names_e0 = {
        {"Gamma.C2.A", -3.0857997869548752234},
        {"Gamma.C2.B", -1.8283034696482216575},
        {"M0.C2.A", -2.0356565020285248835},
        {"M0.C2.B", -2.6544975593746693576},
        {"M1.C2.A", -2.035656502028527548},
        {"M1.C2.B", -2.6544975593746698017},
        {"M2.C2.A", -2.0356565020285266598},
        {"M2.C2.B", -2.6544975593746689135}};
    run_kitaev_gamma_test(K, G, irrep_names_e0);
  }

  {
    double K = -0.70710678118654752441;
    double G = 0.70710678118654752441;
    std::vector<std::pair<std::string, double>> irrep_names_e0 = {
        {"Gamma.C2.A", -3.0140686352095316103},
        {"Gamma.C2.B", -1.7541857925174297872},
        {"M0.C2.A", -1.94808553693512998},
        {"M0.C2.B", -2.5821356185598034472},
        {"M1.C2.A", -1.9480855369351295359},
        {"M1.C2.B", -2.5821356185598021149},
        {"M2.C2.A", -1.948085536935130424},
        {"M2.C2.B", -2.582135618559802559}};
    run_kitaev_gamma_test(K, G, irrep_names_e0);
  }
  Log("done");
}
