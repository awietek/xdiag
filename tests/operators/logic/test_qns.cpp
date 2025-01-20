#include "../../catch.hpp"

#include <xdiag/algebra/matrix.hpp>
#include <xdiag/operators/logic/qns.hpp>

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

} catch (xdiag::Error e) {
  xdiag::error_trace(e);
}
