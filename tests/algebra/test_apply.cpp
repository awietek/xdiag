#include "../catch.hpp"

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/apply.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/algorithms/sparse_diag.hpp>
#include <xdiag/io/read.hpp>
#include <xdiag/operators/logic/order.hpp>
#include <xdiag/operators/logic/symmetrize.hpp>
#include <xdiag/states/fill.hpp>
#include <xdiag/states/random_state.hpp>
#include <xdiag/states/state.hpp>

using namespace xdiag;
using namespace arma;

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

  for (int64_t i = 0; i < nsites; ++i) {
    for (int64_t j = 0; j < nsites; ++j) {
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
        // Log("{} {} {} {} {} {}", i, j, real(m0), real(m1), real(m2), real(m3));
        REQUIRE(isapprox(m0, m2, 1e-6, 1e-6));
        REQUIRE(isapprox(m0, m3, 1e-6, 1e-6));
      }
    }
  }

  int64_t ncols = 10;
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
} catch (Error const &e) {
  error_trace(e);
}
