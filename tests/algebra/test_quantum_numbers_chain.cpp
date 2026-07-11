// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <iostream>

#include <tests/catch.hpp>

#include <xdiag/algebra/symmetrize.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/linalg/sparse_diag.hpp>
#include <xdiag/operators/hc.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/apply.hpp>
#include <xdiag/states/inner.hpp>
#include <xdiag/states/norm.hpp>
#include <xdiag/symmetries/cyclic_group.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;

TEST_CASE("quantum_numbers_chain", "[operators]") try {

  Log("Testing quantum numbers chain");

  for (int nsites = 4; nsites < 10; ++nsites) {
    for (int nup = 1; nup < nsites; ++nup) {

      auto block = Spinhalf(nsites, nup);
      auto H = OpSum();
      for (int s = 0; s < nsites; ++s) {
        H += Op("SdotS", {s, (s + 1) % nsites});
      }
      auto [e0, psi0] = eig0(H, block);

      auto groupi =
          PermutationGroup({Permutation(nsites)}); // identity permutation
      auto irrepi = Representation(groupi);
      auto blocki = Spinhalf(nsites, nup, irrepi);
      auto [e0i, psi0i] = eig0(H, blocki);

      // Determine ground states with momenta
      std::vector<double> e0ks(nsites);
      for (int k = 0; k < nsites; ++k) {
        auto irrep = cyclic_group_irrep(nsites, k);
        auto blockk = Spinhalf(nsites, nup, irrep);
        e0ks[k] = eigval0(H, blockk);
      }
      double e0kmin = *std::min_element(e0ks.begin(), e0ks.end());
      int deg = 0;
      int kmin = 0;
      for (int k = 0; k < nsites; ++k) {
        // Log("k: {}, e0k: {}", k, e0ks[k]);
        if (isapprox(e0ks[k], e0kmin)) {
          ++deg;
          kmin = k;
        }
      }

      REQUIRE(isapprox(e0, e0kmin));
      // Log("nsites: {}, nup: {}, e0: {}, e0k {} (deg: {}, kmin: {})", nsites,
      //     nup, e0, e0kmin, deg, kmin);

      for (int q = 0; q < nsites; ++q) {
        auto kirrep = cyclic_group_irrep(nsites, kmin);
        auto qirrep = cyclic_group_irrep(nsites, q);
        auto Sq = symmetrize(Op("S+", 0), qirrep);
        auto Sq_psi0 = apply(Sq, psi0);
        double nrm0 = norm(Sq_psi0);

        auto Sq_psi0i = apply(Sq, psi0i);
        double nrm0i = norm(Sq_psi0i);
        // Log("q: {}, nrm0: {}, nrm0i: {}", q, nrm0, nrm0i);
        REQUIRE(isapprox(nrm0, nrm0i, 1e-7, 1e-7));

        if (deg == 1) {

          auto blockk = Spinhalf(nsites, nup, kirrep);
          auto [e0k, psi0k] = eig0(H, blockk);
          REQUIRE(isapprox(e0k, e0));

          auto Sq_psi0k = apply(Sq, psi0k);
          auto irr = std::get<Spinhalf>(Sq_psi0k.block()).irreps();
          Representation kqirrep;
          if (isreal(*irr.characters("SitePermutation"))) {
            arma::vec chars =
                (*irr.characters("SitePermutation")).as<arma::vec>();
            kqirrep = Representation(*irr.group("SitePermutation"), chars);
          } else {
            arma::cx_vec chars =
                (*irr.characters("SitePermutation")).as<arma::cx_vec>();
            kqirrep = Representation(*irr.group("SitePermutation"), chars);
          }
          auto kqirrep_exact = cyclic_group_irrep(nsites, (kmin + q) % nsites);
          int nupp1 = *irr.charge("nup");
          REQUIRE(nupp1 == nup + 1);
          REQUIRE(isapprox(kqirrep, kqirrep_exact));

          double nrm0k = norm(Sq_psi0k);
          // Log("q: {}, nrm0: {}, nrm0k: {}", q, nrm0, nrm0k);
          REQUIRE(isapprox(nrm0, nrm0k, 1e-7, 1e-7));
        }
      }
    }
  }

  Log("done");
} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}
