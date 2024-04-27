#include "../catch.hpp"

#include <iostream>

#include "../blocks/spinhalf/testcases_spinhalf.hpp"

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/algorithms/time_evolution/exp_sym_v.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/common.hpp>
#include <xdiag/states/state.hpp>
#include <xdiag/utils/logger.hpp>

TEST_CASE("exp_sym_v", "[algorithms]") try {
  using namespace xdiag;
  for (int64_t N = 2; N <= 10; ++N) {
    Log("Testing exp_sym_v on all-to-all HB model, N={}", N);

    auto block = Spinhalf(N);
    auto psi0 = random_state(block);
    auto bonds = xdiag::testcases::spinhalf::HB_alltoall(N);

    {
      Log("real time");
      // Real time evolution
      double t = 1.2345;
      auto H = matrixC(bonds, block);
      arma::cx_vec psi_ex = expmat(complex(0.0, -1.0 * t) * H) * psi0.vector();
      arma::cx_vec psi = exp_sym_v(bonds, psi0, complex(0, -t)).vectorC();
      // XDIAG_PRINT(norm(psi - psi_ex));
      REQUIRE(norm(psi - psi_ex) < 1e-6);
    }

    {
      Log("imaginary time");
      // Imaginary time evolution
      double t = 1.2345;
      auto H = matrix(bonds, block);
      arma::vec psi_ex = expmat(-t * H) * psi0.vector();
      arma::vec psi = exp_sym_v(bonds, psi0, -t).vector();
      // XDIAG_PRINT(norm(psi - psi_ex));
      REQUIRE(norm(psi - psi_ex) < 1e-6);
    }

    {
      Log("complex time");
      complex t(-1.2345, 3.2145);
      auto H = matrixC(bonds, block);
      arma::cx_vec psi_ex = expmat(t * H) * psi0.vector();
      arma::cx_vec psi = exp_sym_v(bonds, psi0, t).vectorC();
      // XDIAG_PRINT(norm(psi - psi_ex));
      REQUIRE(norm(psi - psi_ex) < 1e-6);
    }
  }

} catch (std::exception const &e) {
  xdiag::traceback(e);
}
