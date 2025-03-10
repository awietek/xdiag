#include <iostream>
#include <string>
#include <vector>
#include <xdiag/all.hpp>

const unsigned N = 36;

auto fl = xdiag::FileToml("triangular.36.60006.J1J2.sublattices.tsl.toml");
toml::table fl2 =
    toml::parse_file("triangular.36.60006.J1J2.sublattices.tsl.toml");
toml::table mat = toml::parse_file("NaYbO2_J.toml");

using namespace xdiag;
void get_OP_Mat(xdiag::OpSum &ops);

int main() try {
  std::vector<std::string> irreps = {
      "Gamma.C2.A", "Gamma.C2.B", "K0.C1.A", "K1.C1.A", "M0.C2.A",
      "M0.C2.B",    "M1.C2.A",    "M1.C2.B", "M2.C2.A", "M2.C2.B",
      "X0.C1.A",    "X1.C1.A",    "X2.C1.A", "Y0.C1.A", "Y1.C1.A",
      "Y2.C1.A",    "Y3.C1.A",    "Y4.C1.A", "Y5.C1.A", "0.C1.A",
      "1.C1.A",     "2.C1.A",     "3.C1.A",  "4.C1.A",  "5.C1.A"};

  arma::mat all_eigen(5, irreps.size(), arma::fill::zeros);

  xdiag::OpSum ops = OpSum();

  get_OP_Mat(ops);
  set_verbosity(2);
  
  for (unsigned i = 0; i < irreps.size(); i++) {
    auto irrep = read_representation(fl, irreps[i], "Symmetries");
    tic();
    auto block = Spinhalf(N, irrep, "3sublattice");
    XDIAG_SHOW(block);
    toc();
    auto res = eigvals_lanczos(ops, block, 1, 1e-12, 1);
    // all_eigen.col(i) = res.eigenvalues.subvec(0, 4);
  }

  // auto save_fl = FileH5("Results/spectrum.h5", "w!");
  // save_fl["spectrum"] = all_eigen;
  return 0;
} catch (xdiag::Error e) {
  xdiag::error_trace(e);
}

void get_OP_Mat(xdiag::OpSum &ops) {
  auto Bonds = fl2["Interactions"];
  unsigned Nb = Bonds.as_array()->size();

  auto sx = arma::cx_mat({{0, std::complex<double>(0.5, 0.0)},
                          {std::complex<double>(0.5, 0.0), 0}});
  auto sy = arma::cx_mat({{0.0, std::complex<double>(0.0, 0.5)},
                          {std::complex<double>(0.0, -0.5), 0.0}});
  auto sz = arma::cx_mat({{std::complex<double>(0.5, 0.0), 0},
                          {0, std::complex<double>(-0.5, 0.0)}});

  arma::cx_mat *matrices[3] = {&sy, &sz, &sx};

  for (unsigned i = 0; i < Nb; ++i) {
    std::string &tag = fl2["Interactions"][i][0].ref<std::string>();
    auto couplingMatrix = mat[tag];

    auto i1 = *fl2["Interactions"][i][2].value<int>();
    auto i2 = *fl2["Interactions"][i][3].value<int>();

    for (unsigned m = 0; m < 3; m++) {
      for (unsigned n = 0; n < 3; n++) {
        arma::cx_mat sasb = arma::kron(*matrices[m], *matrices[n]);
        auto J = *couplingMatrix[m][n].value<double>();
        ops += J * Op("Matrix", {i1, i2}, sasb);
      }
    }
  }
}
