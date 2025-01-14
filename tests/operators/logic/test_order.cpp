#include "../../catch.hpp"

#include <iostream>

#include "../../blocks/electron/testcases_electron.hpp"
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/algorithms/lanczos/eigs_lanczos.hpp>
#include <xdiag/algorithms/sparse_diag.hpp>
#include <xdiag/io/file_toml.hpp>
#include <xdiag/operators/logic/order.hpp>
#include <xdiag/operators/logic/real.hpp>
#include <xdiag/utils/close.hpp>

using namespace xdiag;

TEST_CASE("order", "[operators]") try {
  Log("Testing ordering of operators");
  using namespace arma;

  for (int n_sites = 3; n_sites < 5; ++n_sites) {

    // int n_sites = 6;
    auto ops = testcases::electron::get_linear_chain(n_sites, 1.0, 5.0);
    for (int i = 0; i < n_sites; ++i) {
      ops += "J2" * Op("SdotS", {i, (i + 2) % n_sites});
    }
    ops["J2"] = 0.321;
    ops["T"] = complex(1, 1);

    auto opso = order(ops);
    auto block = Electron(n_sites);
    auto H = matrixC(ops, block);
    auto Ho = matrixC(opso, block);

    REQUIRE(norm(H - Ho) < 1e-12);
  }

  // Kitaev-Gamma test
  {
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/kitaev_gamma/lattice-files/"
                        "honeycomb.8.HeisenbergKitaevGamma.fsl.toml";

    double K = 1.234;
    double G = 3.21;
    auto fl = FileToml(lfile);
    auto ops_read = fl["Interactions"].as<OpSum>();
    auto group = fl["Symmetries"].as<PermutationGroup>();

    cx_mat sx(mat({{0., 0.5}, {0.5, 0.}}), mat({{0., 0.}, {0., 0.}}));
    cx_mat sy(mat({{0., 0.}, {0., 0.}}), mat({{0., -0.5}, {0.5, 0.}}));
    cx_mat sz(mat({{0.5, 0.0}, {0.0, -0.5}}), mat({{0., 0.}, {0., 0.0}}));

    cx_mat sxsx = kron(sx, sx);
    cx_mat sysy = kron(sy, sy);
    cx_mat szsz = kron(sz, sz);
    cx_mat gsx = kron(sy, sz) + kron(sz, sy);
    cx_mat gsy = kron(sx, sz) + kron(sz, sx);
    cx_mat gsz = kron(sx, sy) + kron(sy, sx);

    auto ops = OpSum();
    for (auto [cpl, op] : ops_read) {
      std::string type = op.type();
      auto sites = op.sites();
      if (type == "KITAEVX") {
        ops += K * Op("Matrix", sites, sxsx);
      } else if (type == "KITAEVY") {
        ops += K * Op("Matrix", sites, sysy);
      } else if (type == "KITAEVZ") {
        ops += K * Op("Matrix", sites, szsz);
      } else if (type == "GAMMAX") {
        ops += G * Op("Matrix", sites, gsx);
      } else if (type == "GAMMAY") {
        ops += G * Op("Matrix", sites, gsy);
      } else if (type == "GAMMAZ") {
        ops += G * Op("Matrix", sites, gsz);
      }
    }
    auto opso = order(ops);
    std::vector<std::string> irreps = {"Gamma.C2.A", "Gamma.C2.B", "M0.C2.A",
                                       "M0.C2.B",    "M1.C2.A",    "M1.C2.B",
                                       "M2.C2.A",    "M2.C2.B"};
    for (auto name : irreps) {
      auto irrep = fl[name].as<Representation>();
      auto block = Spinhalf(8, group, irrep);
      cx_mat H = matrixC(ops, block);
      cx_mat Ho = matrixC(opso, block);
      REQUIRE(norm(H - Ho) < 1e-12);
    }
  }

  {
    // test a 3x3 triangular lattice
    Log("tj_symmetric_matrix: tJ 3x3 triangular, order test");
    int64_t n_sites = 9;
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/triangular.9.hop.sublattices.tsl.toml";

    auto fl = FileToml(lfile);
    auto ops = fl["Interactions"].as<OpSum>();
    ops["T"] = 1.0;
    ops["J"] = 0.4;
    auto opso = order(ops);
    auto group = fl["Symmetries"].as<PermutationGroup>();

    // XDIAG_SHOW(ops);
    // XDIAG_SHOW(opso);
    
    std::vector<std::string> rep_name = {
        "Gamma.D3.A1", "Gamma.D3.A2", "Gamma.D3.E", "K0.D3.A1", "K0.D3.A2",
        "K0.D3.E",     "K1.D3.A1",    "K1.D3.A2",   "K1.D3.E",  "Y.C1.A"};

    for (auto name : rep_name) {
      auto irrep = fl[name].as<Representation>();

      for (int64_t nup = 1; nup <= n_sites; ++nup) {
        for (int64_t ndn = 1; ndn <= n_sites; ++ndn) {

          if (nup + ndn > n_sites) {
            continue;
          }
          auto block = tJ(n_sites, nup, ndn, group, irrep);
          cx_mat H = matrixC(ops, block);
          cx_mat Ho = matrixC(opso, block);
          REQUIRE(norm(H - Ho) < 1e-12);
        }
      }
    }
  }

  Log("done");
} catch (xdiag::Error e) {
  xdiag::error_trace(e);
}
