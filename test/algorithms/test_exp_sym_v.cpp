#include "../catch.hpp"

#include <iostream>

#include <hydra/algebra/algebra.h>
#include <hydra/algebra/matrix.h>
#include <hydra/algorithms/time_evolution/exp_sym_v.h>
#include <hydra/blocks/blocks.h>
#include <hydra/common.h>
#include <hydra/states/state.h>
#include <hydra/utils/logger.h>

#include "../blocks/spinhalf/testcases_spinhalf.h"

TEST_CASE("exp_sym_v", "[algorithms]") try {
  using namespace hydra;
  for (int64_t N = 2; N <= 10; ++N) {
    Log("Testing exp_sym_v on all-to-all HB model, N={}", N);

    auto block = Spinhalf(N);
    auto psi0 = random_state(block);
    auto bonds = hydra::testcases::spinhalf::HB_alltoall(N);

    {
      Log("real time");
      // Real time evolution
      double t = 1.2345;
      auto H = matrixC(bonds, block);
      arma::cx_vec psi_ex = expmat(complex(0.0, -1.0 * t) * H) * psi0.vector();
      arma::cx_vec psi = exp_sym_v(bonds, psi0, complex(0, -t)).vectorC();
      // HydraPrint(norm(psi - psi_ex));
      REQUIRE(norm(psi - psi_ex) < 1e-6);
    }

    {
      Log("imaginary time");
      // Imaginary time evolution
      double t = 1.2345;
      auto H = matrix(bonds, block);
      arma::vec psi_ex = expmat(-t * H) * psi0.vector();
      arma::vec psi = exp_sym_v(bonds, psi0, -t).vector();
      // HydraPrint(norm(psi - psi_ex));
      REQUIRE(norm(psi - psi_ex) < 1e-6);
    }

    {
      Log("complex time");
      complex t(-1.2345, 3.2145);
      auto H = matrixC(bonds, block);
      arma::cx_vec psi_ex = expmat(t * H) * psi0.vector();
      arma::cx_vec psi = exp_sym_v(bonds, psi0, t).vectorC();
      // HydraPrint(norm(psi - psi_ex));
      REQUIRE(norm(psi - psi_ex) < 1e-6);
    }
  }

} catch (std::exception const &e) {
  hydra::traceback(e);
}
