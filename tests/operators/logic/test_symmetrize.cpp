// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include <iostream>

#include "../../blocks/electron/testcases_electron.hpp"
#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/apply.hpp>
#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/algorithms/lanczos/eigs_lanczos.hpp>
#include <xdiag/algorithms/sparse_diag.hpp>
#include <xdiag/blocks/electron.hpp>
#include <xdiag/operators/logic/hc.hpp>
#include <xdiag/operators/logic/isapprox.hpp>
#include <xdiag/operators/logic/order.hpp>
#include <xdiag/operators/logic/symmetrize.hpp>
#include <xdiag/states/create_state.hpp>
#include <xdiag/io/read.hpp>
#include <xdiag/algorithms/lanczos/eigvals_lanczos.hpp>

using namespace xdiag;

TEST_CASE("symmetrize", "[operators]") try {
  Log("Testing symmetrize of operator");

  for (int nsites = 3; nsites < 5; ++nsites) {

    // int nsites = 6;
    auto ops = testcases::electron::get_linear_chain(nsites, 1.0, 5.0);
    for (int i = 0; i < nsites; ++i) {
      ops += "J2" * Op("SdotS", {i, (i + 2) % nsites});
    }
    ops["J2"] = 0.321;
    ops["T"] = 0;

    auto irreps = testcases::electron::get_cyclic_group_irreps(nsites);

    for (int nup = 0; nup <= nsites; ++nup) {
      for (int ndn = 0; ndn <= nsites; ++ndn) {
        auto block_nosym = Electron(nsites, nup, ndn);
        if (block_nosym.size() == 0) {
          continue;
        }
        auto o1 = order(ops.plain());
        auto o2 = order(hc(ops).plain());

        // XDIAG_SHOW(o1);
        // XDIAG_SHOW(o2);
        // XDIAG_SHOW(isapprox(o1, o2));

        auto [e0_nosym, v0_nosym] = eig0(ops, block_nosym);

        auto res = eigs_lanczos(ops, block_nosym);
        // XDIAG_SHOW(res.eigenvalues);
        // XDIAG_SHOW(res.criterion);

        {
          auto &v = v0_nosym;
          // Log("v real {}", v.isreal());
          auto e = inner(ops, v);
          // XDIAG_SHOW(norm(v));
          // XDIAG_SHOW(block_nosym);
          // XDIAG_SHOW(e);
          // XDIAG_SHOW(e0_nosym);
          REQUIRE(isapprox(real(e), e0_nosym));
        }
        double e0_sym = 9999.0;
        Representation e0_irrep;
        std::vector<double> e0s;
        for (auto irrep : irreps) {
          auto block_sym = Electron(nsites, nup, ndn, irrep);
          if (block_sym.size() == 0) {
            continue;
          }
          double e0_sector = eigval0(ops, block_sym);

          e0s.push_back(e0_sector);
          if (e0_sector < e0_sym) {
            e0_sym = e0_sector;
            e0_irrep = irrep;
          }
        }

        // Determine degeneracy
        int deg = 0;
        for (auto e0_sector : e0s) {
          if (std::abs(e0_sym - e0_sector) < 1e-6)
            ++deg;
        }
        // Log.out("N: {}, nup:{} ndn: {}, deg: {}, e0_nosym: {}, e0_sym: {}",
        //         nsites, nup, ndn, deg, e0_nosym, e0_sym);
        REQUIRE(isapprox(e0_sym, e0_nosym));

        // Compare correlators only if degeneracy is 1
        // -> g.s. from non-symmetric calculation is unique and symmetric
        if (deg == 1) {
          auto block_sym = Electron(nsites, nup, ndn, e0_irrep);
          // XDIAG_SHOW(ops);

          auto [e0_sym2, v0_sym] = eig0(ops, block_sym);
          REQUIRE(isapprox(e0_sym, e0_nosym));
          {
            auto &v = v0_sym;
            auto Hv = v;
            // Log("UUU");
            apply(ops, v, Hv);
            // Log("VVV");
            auto e = dotC(v, Hv);
            // XDIAG_SHOW(e);
            // XDIAG_SHOW(e0_nosym);
            REQUIRE(isapprox(real(e), e0_nosym));
          }

          // Measure correlators
          for (int i = 1; i < nsites; ++i) {
            OpSum corr_nosym;
            corr_nosym += "J" * Op("SdotS", {0, i});
            corr_nosym["J"] = 1.0;

            auto corr_sym = symmetrize(corr_nosym, e0_irrep.group());
            // Log("XXX");
            auto val_nosym = inner(corr_nosym, v0_nosym);
            // Log("YYY {} {}", corr_sym.isreal(), v0_sym.isreal());
            // XDIAG_SHOW(corr_sym);
            auto val_sym = innerC(corr_sym, v0_sym);
            // Log("ZZZ");
            // Log.out("(0,{}) nosym: {}, sym: {}", i, real(val_nosym),
            //               real(val_sym));
            REQUIRE(isapprox(val_nosym, val_sym, 1e-6, 1e-6));
          }
          // for (int j=0; j<nsites; ++j) {
          //   Ops corr_nosym;
          //   corr_nosym += Op("Exchange", "J", {j, (i+j)%nsites});
          //   Couplings cpls;
          //   cpls["J"] = 1.0;
          //   auto val_nosym2 = Inner(corr_nosym, cpls, block_nosym,
          // v0_nosym);
          //   Log.out("({}, {}) {} {}", j, (i+j)%nsites,
          //   real(val_nosym2),
          //               real(val_nosym));
          // }
        }
        // Log("X");
      }
    }
  }

  OpSum twist;
  twist += complex(0., 1.) * Op("Exchange", {0, 1});
  auto p0 = Permutation({0, 1});
  auto p1 = Permutation({1, 0});
  auto group = PermutationGroup({p0, p1});

  auto twist_sym = symmetrize(twist, group);
  REQUIRE(isapprox(twist_sym, OpSum()));

  auto b = Spinhalf(2);
  auto r = random_state(b);
  auto s = apply(twist_sym, r);
  REQUIRE(isapprox(dot(r, s), 0.));

  {
    Log("Leo's bug report Electron #85");
    std::string latticeInput =
        XDIAG_DIRECTORY "/misc/data/hubbard.cluster.8.toml";
    auto lfile = FileToml(latticeInput);

    // Creating the Hilbert
    int N = 8;
    int nup = 4;
    int ndn = 4;

    // Creating the Hamiltonian
    auto ham = read_opsum(lfile, "Interactions");

    double t = 0.436;
    double tp = -0.07;
    double tpp = 0.05;
    double U = 3.12;

    ham["Tx"] = t;
    ham["Ty"] = t;
    ham["Tp+"] = tp;
    ham["Tp-"] = tp;
    ham["Tppx"] = tpp;
    ham["Tppy"] = tpp;
    ham["U"] = U;

    // Creating the irreps
    std::vector<Representation> irreps;
    for (unsigned int i = 0; i < 8; ++i) {
      auto irrep = read_representation(lfile, fmt::format("Irrep{}", i));
      irreps.push_back(irrep);
    }

    // Ground state is at Gamma
    auto irrep = irreps[0];
    auto block = Electron(N, nup, ndn, irrep);
    auto block_nosym = Electron(N, nup, ndn);

    auto [e0, gs] = eig0(ham, block);
    auto [e0n, gsn] = eig0(ham, block_nosym);

    Log("Found ground state energy : {:10.6f} {}", e0, e0n);
    gs.make_complex();
    gsn.make_complex();

    for (unsigned int i = 0; i < 8; ++i) {

      std::cout << "Irrep " << i << std::endl;
      auto S_q = symmetrize(Op("Sz", 0), irreps[i]);
      std::cout << "symmetrized.." << std::endl;
      auto Av = apply(S_q, gs);
      auto Avn = apply(S_q, gsn);

      complex e = innerC(ham, Av);
      complex en = innerC(ham, Avn);
      // Log("check: {} {}", e, en);
      REQUIRE(isapprox(e, en, 1e-6, 1e-6));
      std::cout << "applied.." << std::endl;
      auto b = Av.block();
      auto nrm = norm(Av);
      Av /= nrm;

      auto res = eigvals_lanczos_inplace(ham, Av);
    }
  }

  Log("done");
} catch (xdiag::Error e) {
  xdiag::error_trace(e);
}
