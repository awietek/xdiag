// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../catch.hpp"

#include <xdiag/algebra/matrix.hpp>
#include <xdiag/blocks/electron.hpp>
#include <xdiag/common.hpp>

using namespace xdiag;
using namespace arma;

TEST_CASE("algebra_matrix", "[algebra]") try {
  Log("Test algebra matrix");

  for (int64_t nsites = 2; nsites < 5; ++nsites) {
    for (int64_t nup = 0; nup <= nsites; ++nup) {
      for (int64_t ndn = 0; ndn <= nsites; ++ndn) {
        auto block = Electron(nsites, nup, ndn);

        for (int i = 0; i < nsites; ++i) {
          for (int j = 0; j < nsites; ++j) {
            if (i != j) {
              auto mi = matrix(Op("Ntot", i), block);
              auto mj = matrix(Op("Ntot", j), block);
              auto mij = matrix(Op("NtotNtot", {i, j}), block);
              REQUIRE(norm(mij - mi * mj) < 1e-12);
            }

            if ((nup > 0) && (nup < nsites)) {
              auto block0 = Electron(nsites, nup, ndn);
              auto blockp = Electron(nsites, nup + 1, ndn);
              auto blockm = Electron(nsites, nup - 1, ndn);

              auto cdagupir = matrix(Op("Cdagup", i), block0);
              auto cupjr = matrix(Op("Cup", j), blockp);

              auto cupjl = matrix(Op("Cup", j), block0);
              auto cdagupil = matrix(Op("Cdagup", i), blockm);

              arma::mat anti_comm = cdagupil * cupjl + cupjr * cdagupir;
              int64_t D = block0.size();
              mat id(D, D, fill::eye);

              if (i == j) {
                REQUIRE(norm(anti_comm - id) < 1e-12);
              } else {
                REQUIRE(norm(anti_comm) < 1e-12);
              }
            }
          }
        }
      }
    }
  }
} catch (Error const &e) {
  error_trace(e);
}
