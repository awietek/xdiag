// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <mpi.h>

#include <tests/blocks/electron/testcases_electron.hpp>
#include <tests/blocks/spinhalf/testcases_spinhalf.hpp>
#include <tests/blocks/tj/testcases_tj.hpp>
#include <tests/catch.hpp>

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/blocks/distributed/tj_distributed.hpp>
#include <xdiag/blocks/tj.hpp>
#include <xdiag/linalg/sparse_diag.hpp>
#include <xdiag/matrices/apply.hpp>
#include <xdiag/matrices/matrix.hpp>
#include <xdiag/states/apply.hpp>
#include <xdiag/states/create_state.hpp>
#include <xdiag/states/fill.hpp>
#include <xdiag/states/norm.hpp>
#include <xdiag/states/dot.hpp>
#include <xdiag/states/inner.hpp>
#include <xdiag/utils/logger.hpp>

TEST_CASE("electron_distributed_raise_lower", "[electron_distributed]") try {
  // using block_t = Electron;
  using namespace xdiag;
  using namespace testcases::electron;
  using block_t = ElectronDistributed;

  // Test normal ordering
  Log("test ElectronDistributed normal ordering");
  for (int nsites = 2; nsites < 6; ++nsites) {
    auto block0 = ElectronDistributed(nsites, 0, 0);
    auto psi0 = product_state(block0, std::vector<int64_t>(nsites, 0));

    for (int nup = 0; nup <= nsites; ++nup) {
      for (int ndn = 0; ndn <= nsites; ++ndn) {
        auto block = ElectronDistributed(nsites, nup, ndn);
        auto blockserial = Electron(nsites, nup, ndn);

        for (auto pstate : blockserial) {
          std::vector<int> up_positions;
          std::vector<int> dn_positions;
          for (int i = 0; i < nsites; ++i) {
            if ((pstate[i] == 1) || (pstate[i] == 3)) {
              up_positions.push_back(i);
            }
            if ((pstate[i] == 2) || (pstate[i] == 3)) {
              dn_positions.push_back(i);
            }
          }

          // Create state from product state
          auto psi = State(block);
          fill(psi, pstate);

          auto psi2 = psi0;

          std::reverse(dn_positions.begin(), dn_positions.end());
          for (int i : dn_positions) {
            psi2 = apply(Op("Cdagdn", i), psi2);
          }

          std::reverse(up_positions.begin(), up_positions.end());
          for (int i : up_positions) {
            psi2 = apply(Op("Cdagup", i), psi2);
          }

          REQUIRE(norm(psi2 - psi) < 1e-12);
        }
      }
    }
  }

  using namespace arma;

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
