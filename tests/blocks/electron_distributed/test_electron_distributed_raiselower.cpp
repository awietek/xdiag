// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"
#include "../electron/testcases_electron.hpp"
#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/apply.hpp>
#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/states/create_state.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/xdiag_show.hpp>

using namespace xdiag;
using namespace arma;

TEST_CASE("electron_distributed_raise_lower", "[electron_distributed]") try {
  using namespace testcases::electron;
  using block_t = ElectronDistributed;
  // using block_t = Electron;

  std::vector<std::string> op_strs = {"Cdagup", "Cdagdn", "Cup", "Cdn"};
  // std::vector<std::string> op_strs = {"Cup", "Cdn"};

  for (int nsites = 2; nsites < 6; ++nsites) {
    Log("testing electron anticommutation relations: N={}", nsites);

    for (int nup = 0; nup <= nsites; ++nup) {
      for (int ndn = 0; ndn <= nsites; ++ndn) {

        auto block = block_t(nsites, nup, ndn);

        for (int i = 0; i < nsites; ++i) {
          for (int j = 0; j < nsites; ++j) {

            for (auto op_i_str : op_strs) {
              for (auto op_j_str : op_strs) {

                if (!valid_nup_ndn(op_i_str, op_j_str, nup, ndn, nsites)) {
                  continue;
                }
                if (!valid_nup_ndn(op_i_str, nup, ndn, nsites)) {
                  continue;
                }
                if (!valid_nup_ndn(op_j_str, nup, ndn, nsites)) {
                  continue;
                }
                // Log("nsites: {}, nup: {}, ndn: {}, i: {}, j: {}, op_i: {}, "
                //     "op_j: {} ",
                //     nsites, nup, ndn, i, j, op_i_str, op_j_str);

                auto op_i = Op(op_i_str, i);
                auto op_j = Op(op_j_str, j);
                auto r = random_state(block);
                auto anti_comm =
                    apply(op_i, apply(op_j, r)) + apply(op_j, apply(op_i, r));

                // check the fermionic commutation relations
                if (i != j) {
                  REQUIRE(norm(anti_comm) < 1e-12);
                } else {
                  if (((op_i_str == "Cdagup") && (op_j_str == "Cup")) ||
                      ((op_i_str == "Cup") && (op_j_str == "Cdagup")) ||
                      ((op_i_str == "Cdagdn") && (op_j_str == "Cdn")) ||
                      ((op_i_str == "Cdn") && (op_j_str == "Cdagdn"))) {
                    REQUIRE(norm(anti_comm - r) < 1e-12);
                  } else {
                    REQUIRE(norm(anti_comm) < 1e-12);
                  }
                }

              } // loop op_i_str
            } // loop op_j_str

            // Check hopping
            if ((nup > 0) && (nup < nsites)) {
              auto r = random_state(block);
              auto a = apply(Op("Hopup", {i, j}), r);
              auto b = -(apply(Op("Cdagup", i), apply(Op("Cup", j), r)) +
                         apply(Op("Cdagup", j), apply(Op("Cup", i), r)));

              REQUIRE(isapprox(a, b));
            }

          } // loop i
        } // loop j

      } // loop nup
    } // loop ndn

  } // loop nsites
} catch (xdiag::Error const &e) {
  error_trace(e);
}
