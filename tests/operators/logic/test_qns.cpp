// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include "../../blocks/electron/testcases_electron.hpp"
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/io/read.hpp>
#include <xdiag/operators/logic/order.hpp>
#include <xdiag/operators/logic/qns.hpp>
#include <xdiag/operators/logic/symmetrize.hpp>
#include <xdiag/utils/xdiag_show.hpp>

using namespace xdiag;
using namespace arma;

TEST_CASE("qns", "[operators]") try {
  Log("Testing quantum number computations of Ops and OpSums");
  REQUIRE(*nup(Op("S+", 0)) == 1);
  REQUIRE(*ndn(Op("S+", 0)) == -1);
  REQUIRE(*nup(Op("S-", 0)) == -1);
  REQUIRE(*ndn(Op("S-", 0)) == 1);

  REQUIRE(*nup(Op("Cdagup", 0)) == 1);
  REQUIRE(*ndn(Op("Cdagup", 0)) == 0);
  REQUIRE(*nup(Op("Cup", 0)) == -1);
  REQUIRE(*ndn(Op("Cup", 0)) == 0);

  REQUIRE(*nup(Op("Cdagdn", 0)) == 0);
  REQUIRE(*ndn(Op("Cdagdn", 0)) == 1);
  REQUIRE(*nup(Op("Cdn", 0)) == 0);
  REQUIRE(*ndn(Op("Cdn", 0)) == -1);

  REQUIRE(*nup(Op("SdotS", {0, 1})) == 0);
  REQUIRE(*ndn(Op("SdotS", {0, 1})) == 0);
  REQUIRE(*nup(Op("Exchange", {0, 1})) == 0);
  REQUIRE(*ndn(Op("Exchange", {0, 1})) == 0);
  REQUIRE(*nup(Op("Sz", 0)) == 0);
  REQUIRE(*ndn(Op("Sz", 0)) == 0);
  REQUIRE(*nup(Op("SzSz", {0, 1})) == 0);
  REQUIRE(*ndn(Op("SzSz", {0, 1})) == 0);
  REQUIRE(*nup(Op("ScalarChirality", {0, 1, 2})) == 0);
  REQUIRE(*ndn(Op("ScalarChirality", {0, 1, 2})) == 0);
  REQUIRE(*nup(Op("Hop", {0, 1})) == 0);
  REQUIRE(*ndn(Op("Hop", {0, 1})) == 0);
  REQUIRE(*nup(Op("Hopup", {0, 1})) == 0);
  REQUIRE(*ndn(Op("Hopup", {0, 1})) == 0);
  REQUIRE(*nup(Op("Hopdn", {0, 1})) == 0);
  REQUIRE(*ndn(Op("Hopdn", {0, 1})) == 0);
  REQUIRE(*nup(Op("tJSdotS", {0, 1})) == 0);
  REQUIRE(*ndn(Op("tJSdotS", {0, 1})) == 0);
  REQUIRE(*nup(Op("tJSzSz", {0, 1})) == 0);
  REQUIRE(*ndn(Op("tJSzSz", {0, 1})) == 0);
  REQUIRE(*nup(Op("Nup", 0)) == 0);
  REQUIRE(*ndn(Op("Nup", 0)) == 0);
  REQUIRE(*nup(Op("Ntot", 0)) == 0);
  REQUIRE(*ndn(Op("Ntot", 0)) == 0);
  REQUIRE(*nup(Op("Ndn", 0)) == 0);
  REQUIRE(*ndn(Op("Ndn", 0)) == 0);
  REQUIRE(*nup(Op("HubbardU")) == 0);
  REQUIRE(*ndn(Op("HubbardU")) == 0);

  mat sx({{0., 0.5}, {0.5, 0.}});
  cx_mat sy(mat({{0., 0.}, {0., 0.}}), mat({{0., -0.5}, {0.5, 0.}}));
  mat sz({{0.5, 0.0}, {0.0, -0.5}});
  mat sm({{0., 1.}, {0., 0.}});
  mat sp({{0., 0.}, {1., 0.}});

  REQUIRE(nup(Op("Matrix", 0, sx)) == std::nullopt);
  REQUIRE(nup(Op("Matrix", 0, sy)) == std::nullopt);
  REQUIRE(*nup(Op("Matrix", 0, sz)) == 0);
  REQUIRE(*nup(Op("Matrix", 0, sp)) == 1);
  REQUIRE(*nup(Op("Matrix", 0, sm)) == -1);

  REQUIRE(*nup(Op("Matrix", {0, 1}, mat(kron(sp, sp)))) == 2);
  REQUIRE(*nup(Op("Matrix", {0, 1}, mat(kron(sp, sz)))) == 1);
  REQUIRE(*nup(Op("Matrix", {0, 1}, mat(kron(sp, sm)))) == 0);

  REQUIRE(*nup(Op("Matrix", {0, 1}, mat(kron(sz, sp)))) == 1);
  REQUIRE(*nup(Op("Matrix", {0, 1}, mat(kron(sz, sz)))) == 0);
  REQUIRE(*nup(Op("Matrix", {0, 1}, mat(kron(sz, sm)))) == -1);

  REQUIRE(*nup(Op("Matrix", {0, 1}, mat(kron(sm, sp)))) == 0);
  REQUIRE(*nup(Op("Matrix", {0, 1}, mat(kron(sm, sz)))) == -1);
  REQUIRE(*nup(Op("Matrix", {0, 1}, mat(kron(sm, sm)))) == -2);

  for (int nsites = 3; nsites < 8; ++nsites) {
    auto irreps = testcases::electron::get_cyclic_group_irreps(nsites);
    for (auto const &irrep : irreps) {
      {
        auto ops = symmetrize(Op("Sz", 0), irrep);
        auto irrep2 = representation(ops, irrep.group());
        REQUIRE(irrep2);
        REQUIRE(isapprox(irrep, *irrep2));
      }

      {
        auto opss = OpSum();
        opss += Op("SdotS", {0, 1});
        opss += Op("SdotS", {1, 2});
        auto ops = symmetrize(Op("Sz", 0), irrep);
        auto irrep2 = representation(ops, irrep.group());
        REQUIRE(irrep2);
        REQUIRE(isapprox(irrep, *irrep2));
      }
    }
  }

  std::string lfile = XDIAG_DIRECTORY
      "/misc/data/triangular.9.tup.phi.tdn.nphi.sublattices.tsl.toml";

  auto fl = FileToml(lfile);
  auto ops = read_opsum(fl, "Interactions");
  double eta = 0.1;
  ops["TPHI"] = complex(cos(eta * M_PI), sin(eta * M_PI));
  ops["JPHI"] = complex(cos(2 * eta * M_PI), sin(2 * eta * M_PI));
  auto irrep = read_representation(fl, "Gamma.D3.A1");
  auto irrep2 = representation(ops, irrep.group());
  REQUIRE(irrep2);
  REQUIRE(isapprox(*irrep2, irrep));

  // Testing combined qns of matrices (issue #17)
  {
    std::complex<double> imaginary_unit{0, 1};
    auto sx = arma::mat({{0, 1}, {1, 0}});
    auto sy = arma::cx_mat({{0, -imaginary_unit}, {imaginary_unit, 0}});
    auto sz = arma::mat({{1.0, 0.0}, {0.0, -1.0}});

    arma::mat sxsx = arma::kron(sx, sx);
    arma::cx_mat sysy = arma::kron(sy, sy);
    arma::mat szsz = arma::kron(sz, sz);

    auto SdotS = OpSum();
    SdotS += Op("Matrix", {0, 1}, sxsx);
    SdotS += Op("Matrix", {0, 1}, sysy);
    SdotS += Op("Matrix", {0, 1}, szsz);

    REQUIRE(matrixC(SdotS, Spinhalf(2, 1)).n_cols == 2);
  }

} catch (xdiag::Error e) {
  xdiag::error_trace(e);
}
