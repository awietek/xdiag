// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../catch.hpp"

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/apply.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/algorithms/sparse_diag.hpp>
#include <xdiag/io/read.hpp>
#include <xdiag/operators/logic/order.hpp>
#include <xdiag/operators/logic/real.hpp>
#include <xdiag/operators/logic/symmetrize.hpp>
#include <xdiag/states/create_state.hpp>
#include <xdiag/states/fill.hpp>
#include <xdiag/states/random_state.hpp>
#include <xdiag/states/state.hpp>

using namespace xdiag;
using namespace arma;

void test_apply_matrix(OpSum ops, Block block, bool real, int64_t ncols) try {
  int64_t seed = 42;
  bool normalized = true;

  auto r = random_state(block, real, ncols, seed, normalized);

  auto Ar = apply(ops, r);

  if (real && isreal(ops)) {
    auto rAr = matrix_dot(r, Ar);
    for (int i = 0; i < ncols; ++i) {
      for (int j = 0; j < ncols; ++j) {
        double rAr_ij = dot(r.col(i), apply(ops, r.col(j)));
        REQUIRE(isapprox(rAr_ij, rAr(i, j)));
      }
    }
  } else {

    auto rAr = matrix_dotC(r, Ar);
    for (int i = 0; i < ncols; ++i) {
      for (int j = 0; j < ncols; ++j) {
        complex rAr_ij = dotC(r.col(i), apply(ops, r.col(j)));
        REQUIRE(isapprox(rAr_ij, rAr(i, j)));
      }
    }
  }
}
XDIAG_CATCH

TEST_CASE("algebra_apply", "[algebra]") try {
  Log("Test algebra apply");
  std::string lfile = XDIAG_DIRECTORY
      "/misc/data/triangular.9.tup.phi.tdn.nphi.sublattices.tsl.toml";
  auto fl = FileToml(lfile);
  auto ops = read_opsum(fl, "Interactions");
  ops += "U" * Op("HubbardU");
  ops["TPHI"] = complex(0.5, 0.1);
  ops["JPHI"] = 0.;
  ops["U"] = 5.;
  // std::vector<std::string> irreps = {
  //     "Gamma.D3.A1", "Gamma.D3.A2", "Gamma.D3.E", "K0.D3.A1", "K0.D3.A2",
  //     "K0.D3.E",     "K1.D3.A1",    "K1.D3.A2",   "K1.D3.E",  "Y.C1.A"};

  int64_t nsites = 9;
  int64_t nup = 4;
  int64_t ndn = 4;

  auto block = Electron(nsites, nup, ndn);
  auto [e0, gs] = eig0(ops, block);

  auto irrep = read_representation(fl, "Gamma.D3.A2");
  auto blocksym = Electron(nsites, nup, ndn, irrep);
  auto [e0sym, gssym] = eig0(ops, blocksym);
  REQUIRE(isapprox(e0, e0sym, 1e-8, 1e-8));
  // Log("e0 {} {}", e0, e0sym);

  for (int64_t i = 0; i < nsites / 2; ++i) {
    for (int64_t j = 0; j < nsites / 2; ++j) {
      auto nigs = apply(Op("Ntot", i), gs);
      auto njgs = apply(Op("Ntot", j), gs);
      complex m0 = dotC(nigs, njgs);
      if (i != j) {
        complex m1 = innerC(Op("NtotNtot", {i, j}), gs);
        REQUIRE(isapprox(m0, m1));
        auto nijs = symmetrize(Op("NtotNtot", {i, j}), irrep.group());
        // XDIAG_SHOW(order(nijs));
        complex m2 = innerC(nijs, gssym);
        auto nijsgs = apply(nijs, gssym);
        complex m3 = dotC(gssym, nijsgs);
        // Log("{} {} {} {} {} {}", i, j, real(m0), real(m1), real(m2),
        // real(m3));
        REQUIRE(isapprox(m0, m2, 1e-6, 1e-6));
        REQUIRE(isapprox(m0, m3, 1e-6, 1e-6));
      }

      auto niu_gs = apply(Op("Nup", i), gs);
      auto nid_gs = apply(Op("Ndn", i), gs);
      auto nju_gs = apply(Op("Nup", j), gs);
      auto njd_gs = apply(Op("Ndn", j), gs);

      complex m_iujd = dotC(niu_gs, njd_gs);
      complex m_juid = dotC(nju_gs, nid_gs);
      complex m_uu = dotC(niu_gs, nju_gs);
      complex m_dd = dotC(nid_gs, njd_gs);

      if (i != j) {
        complex m1 = innerC(Op("NupNdn", {i, j}), gs);
        REQUIRE(isapprox(m_iujd, m1));
        m1 = innerC(Op("NupNdn", {j, i}), gs);
        REQUIRE(isapprox(m_juid, m1));
        m1 = innerC(Op("NupNup", {i, j}), gs);
        REQUIRE(isapprox(m_uu, m1));
        m1 = innerC(Op("NdnNdn", {i, j}), gs);
        REQUIRE(isapprox(m_dd, m1));

        auto niujd_gs_s = symmetrize(Op("NupNdn", {i, j}), irrep.group());
        complex m2 = innerC(niujd_gs_s, gssym);
        auto niujd_s_gs = apply(niujd_gs_s, gssym);
        complex m3 = dotC(gssym, niujd_s_gs);
        REQUIRE(isapprox(m_iujd, m2, 1e-6, 1e-6));
        REQUIRE(isapprox(m_iujd, m3, 1e-6, 1e-6));

        auto njuid_gs_s = symmetrize(Op("NupNdn", {j, i}), irrep.group());
        m2 = innerC(njuid_gs_s, gssym);
        auto njuid_s_gs = apply(njuid_gs_s, gssym);
        m3 = dotC(gssym, njuid_s_gs);
        REQUIRE(isapprox(m_juid, m2, 1e-6, 1e-6));
        REQUIRE(isapprox(m_juid, m3, 1e-6, 1e-6));

        auto nuu_gs_s = symmetrize(Op("NupNup", {i, j}), irrep.group());
        m2 = innerC(nuu_gs_s, gssym);
        auto nuu_s_gs = apply(nuu_gs_s, gssym);
        m3 = dotC(gssym, nuu_s_gs);
        REQUIRE(isapprox(m_uu, m2, 1e-6, 1e-6));
        REQUIRE(isapprox(m_uu, m3, 1e-6, 1e-6));

        auto ndd_gs_s = symmetrize(Op("NdnNdn", {i, j}), irrep.group());
        m2 = innerC(ndd_gs_s, gssym);
        auto ndd_s_gs = apply(ndd_gs_s, gssym);
        m3 = dotC(gssym, ndd_s_gs);
        REQUIRE(isapprox(m_dd, m2, 1e-6, 1e-6));
        REQUIRE(isapprox(m_dd, m3, 1e-6, 1e-6));
      }
    }
  }

  int64_t ncols = 3;
  auto v = State(block, false, ncols);
  auto vs = State(blocksym, false, ncols);
  for (int64_t m = 0; m < ncols; m++) {
    auto rstate = RandomState(m, false);
    xdiag::fill(v, rstate, m);
    xdiag::fill(vs, rstate, m);
  }

  auto Hv = apply(ops, v);
  auto Hvs = apply(ops, vs);

  auto hmat = xdiag::matrixC(ops, block);
  auto hmats = xdiag::matrixC(ops, blocksym);

  for (int64_t m = 0; m < ncols; ++m) {
    auto evec = v.vectorC(m, false);

    auto n0 = as_scalar(evec.t() * hmat * evec);
    auto n1 = dotC(v.col(m), Hv.col(m));
    REQUIRE(isapprox(n0, n1, 1e-6, 1e-6));

    auto evecsym = vs.vectorC(m, false);

    auto n0s = as_scalar(evecsym.t() * hmats * evecsym);
    auto n1s = dotC(vs.col(m), Hvs.col(m));
    REQUIRE(isapprox(n0s, n1s, 1e-6, 1e-6));
  }

  {
    OpSum ops;
    ops += Op("Hop", {0, 1});
    Electron block(2);
    test_apply_matrix(ops, block, true, 3);
    test_apply_matrix(ops, block, true, 4);
    test_apply_matrix(ops, block, false, 3);
    test_apply_matrix(ops, block, false, 4);
  }

  {
    OpSum ops;
    ops += Op("Cdagup", 1);
    Electron block(2);
    test_apply_matrix(ops, block, true, 3);
    test_apply_matrix(ops, block, true, 4);
    test_apply_matrix(ops, block, false, 3);
    test_apply_matrix(ops, block, false, 4);
  }
} catch (Error const &e) {
  error_trace(e);
}
