#include "../../catch.hpp"
#include "../tj/testcases_tj.hpp"
#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/apply.hpp>
#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/states/create_state.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;
using namespace arma;

TEST_CASE("tj_distributed_raise_lower", "[tj_distributed]") try {
  using namespace testcases::tj;
  using block_t = tJDistributed;
  // using block_t = tJ;

  std::vector<std::string> op_strs = {"Cdagup", "Cdagdn", "Cup", "Cdn"};
  // std::vector<std::string> op_strs = {"Cup", "Cdn"};

  for (int nsites = 2; nsites < 6; ++nsites) {
    Log("testing tj anticommutation relations: N={}", nsites);

    for (int nup = 0; nup <= nsites; ++nup) {
      for (int ndn = 0; ndn <= nsites - nup; ++ndn) {

        auto block = block_t(nsites, nup, ndn);

        for (int i = 0; i < nsites; ++i) {
          for (int j = 0; j < nsites; ++j) {

            for (auto op_i_str : op_strs) {
              for (auto op_j_str : op_strs) {
                // Log("nsites: {}, nup: {}, ndn: {}, i: {}, j: {}, op_i: {}, "
                //     "op_j: {} ",
                //     nsites, nup, ndn, i, j, op_i_str, op_j_str);

                if (!valid_nup_ndn(op_i_str, op_j_str, nup, ndn, nsites)) {
                  continue;
                }
                if (!valid_nup_ndn(op_i_str, nup, ndn, nsites)) {
                  continue;
                }
                if (!valid_nup_ndn(op_j_str, nup, ndn, nsites)) {
                  continue;
                }

                auto op_i = Op(op_i_str, i);
                auto op_j = Op(op_j_str, j);

                auto r = random_state(block);
                auto anti_comm =
                    apply(op_i, apply(op_j, r)) + apply(op_j, apply(op_i, r));

                // check the non-fermionic commutation relations of the t-J
                // model
                // see arxiv.org/abs/0706,4236 (tJ model then and now ... Jozef
                // Spalek)
                if (i != j) {
                  REQUIRE(norm(anti_comm) < 1e-12);
                } else {
                  if (((op_i_str == "Cdagup") && (op_j_str == "Cdagup")) ||
                      ((op_i_str == "Cdagup") && (op_j_str == "Cdagdn")) ||
                      ((op_i_str == "Cdagdn") && (op_j_str == "Cdagup")) ||
                      ((op_i_str == "Cdagdn") && (op_j_str == "Cdagdn")) ||
                      ((op_i_str == "Cup") && (op_j_str == "Cup")) ||
                      ((op_i_str == "Cup") && (op_j_str == "Cdn")) ||
                      ((op_i_str == "Cdn") && (op_j_str == "Cup")) ||
                      ((op_i_str == "Cdn") && (op_j_str == "Cdn"))) {
                    REQUIRE(norm(anti_comm) < 1e-12);
                  } else if (((op_i_str == "Cdagup") && (op_j_str == "Cup")) ||
                             ((op_i_str == "Cup") && (op_j_str == "Cdagup"))) {
                    auto ndn_op = Op("Ndn", i);
                    REQUIRE(norm(anti_comm + apply(ndn_op, r) - r) < 1e-12);
                  } else if (((op_i_str == "Cdagdn") && (op_j_str == "Cdn")) ||
                             ((op_i_str == "Cdn") && (op_j_str == "Cdagdn"))) {
                    auto nup_op = Op("Nup", i);
                    REQUIRE(norm(anti_comm + apply(nup_op, r) - r) < 1e-12);
                  } else if ((op_i_str == "Cdagup") && (op_j_str == "Cdn")) {
                    REQUIRE(norm(anti_comm - apply(Op(op_i_str, i),
                                                   apply(Op(op_j_str, j), r))) <
                            1e-12);
                  } else if ((op_i_str == "Cdn") && (op_j_str == "Cdagup")) {
                    REQUIRE(norm(anti_comm - apply(Op(op_j_str, j),
                                                   apply(Op(op_i_str, i), r))) <
                            1e-12);
                  } else if ((op_i_str == "Cdagdn") && (op_j_str == "Cup")) {
                    REQUIRE(norm(anti_comm - apply(Op(op_i_str, i),
                                                   apply(Op(op_j_str, j), r))) <
                            1e-12);

                  } else if ((op_i_str == "Cup") && (op_j_str == "Cdagdn")) {
                    REQUIRE(norm(anti_comm - apply(Op(op_j_str, j),
                                                   apply(Op(op_i_str, i), r))) <
                            1e-12);
                  } else {
                    Log("unchecked"); // shouldn't happen
                  }
                }

              } // loop op_i_str
            } // loop op_j_str

          } // loop i
        } // loop j

      } // loop nup
    } // loop ndn

  } // loop nsites
} catch (xdiag::Error const &e) {
  error_trace(e);
}
