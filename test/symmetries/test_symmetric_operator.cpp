#include "../catch.hpp"

#include <iostream>

#include "../blocks/electron/testcases_electron.h"
#include <hydra/all.h>

using namespace hydra;

TEST_CASE("symmetric_operator", "[symmetries]") {
  lila::Log("Testing symmetric_operator");

  for (int n_sites = 3; n_sites < 5; ++n_sites) {

    // int n_sites = 6;
    auto [bondlist, couplings] =
        testcases::electron::get_linear_chain(n_sites, 1.0, 5.0);
    for (int i = 0; i < n_sites; ++i) {
      bondlist << Bond("HB", "J2", {i, (i + 2) % n_sites});
    }
    couplings["J2"] = 0.321;

    auto [space_group, irreps] =
        testcases::electron::get_cyclic_group_irreps(n_sites);

    for (int nup = 0; nup <= n_sites; ++nup) {
      for (int ndn = 0; ndn <= n_sites; ++ndn) {

        auto block_nosym = Electron(n_sites, nup, ndn);
        if (block_nosym.size() == 0)
          continue;

        auto [e0_nosym, v0_nosym] =
            GroundstateCplx(bondlist, couplings, block_nosym);

        {
          auto &v = v0_nosym;
          auto Hv = v;
          Apply(bondlist, couplings, v, Hv);
          auto e = Dot(v, Hv);
          REQUIRE(lila::close(lila::real(e), e0_nosym));
        }
        double e0_sym = 9999.0;
        Representation e0_irrep;
        std::vector<double> e0s;
        for (auto irrep : irreps) {
          auto block_sym = Electron(n_sites, nup, ndn, space_group, irrep);
          if (block_sym.size() == 0)
            continue;
          double e0_sector = E0Cplx(bondlist, couplings, block_sym);
          e0s.push_back(e0_sector);
          // lila::Log.out("sector {}", e0_sector);
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
        // lila::Log.out(
        //     "N: {}, nup:{} ndn: {}, deg: {}, e0_nosym: {}, e0_sym: {}",
        //     n_sites, nup, ndn, deg, e0_nosym, e0_sym);
        REQUIRE(lila::close(e0_sym, e0_nosym));

        // Compare correlators only if degeneracy is 1
        // -> g.s. from non-symmetric calculation is unique and symmetric
        if (deg == 1) {
          auto block_sym = Electron(n_sites, nup, ndn, space_group, e0_irrep);
          auto [e0_sym2, v0_sym] =
              GroundstateCplx(bondlist, couplings, block_sym);
          REQUIRE(lila::close(e0_sym, e0_nosym));
          {
            auto &v = v0_sym;
            auto Hv = v;
            Apply(bondlist, couplings, v, Hv);
            auto e = Dot(v, Hv);
            REQUIRE(lila::close(lila::real(e), e0_nosym));
          }
          // Measure correlators
          for (int i = 1; i < n_sites; ++i) {
            BondList corr_nosym;
            corr_nosym << Bond("HB", "J", {0, i});
            Couplings cpls;
            cpls["J"] = 1.0;
            auto [corr_sym, cpls_sym] =
                SymmetricOperator(corr_nosym, cpls, space_group);

            auto val_nosym = Inner(corr_nosym, cpls, v0_nosym);
            auto val_sym = Inner(corr_sym, cpls_sym, v0_sym);
            // lila::Log.out("(0,{}) nosym: {}, sym: {}", i,
            // lila::real(val_nosym),
            //               lila::real(val_sym));
            REQUIRE(lila::close(val_nosym, val_sym, 1e-6, 1e-6));
          }
          // for (int j=0; j<n_sites; ++j) {
          //   BondList corr_nosym;
          //   corr_nosym << Bond("EXCHANGE", "J", {j, (i+j)%n_sites});
          //   Couplings cpls;
          //   cpls["J"] = 1.0;
          //   auto val_nosym2 = Inner(corr_nosym, cpls, block_nosym, v0_nosym);
          //   lila::Log.out("({}, {}) {} {}", j, (i+j)%n_sites,
          //   lila::real(val_nosym2),
          //               lila::real(val_nosym));
          // }
        }
      }
    }
  }
  lila::Log("done");
}
