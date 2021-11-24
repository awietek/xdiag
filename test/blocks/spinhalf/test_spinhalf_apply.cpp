#include "../../catch.hpp"

#include <iostream>

#include "testcases_spinhalf.h"
#include <hydra/all.h>
#include <lila/all.h>

using namespace hydra;

void test_apply(BondList bonds, Couplings couplings){
  int N = bonds.n_sites();
  for (int nup = 0; nup <= N; ++nup) {
    auto block = Spinhalf<uint32_t>(N, nup);
    auto H = MatrixReal(bonds, couplings, block, block);
    REQUIRE(lila::close(H, lila::Herm(H)));
    
    auto v = lila::Random<double>(block.size());
    auto w1 = lila::Mult(H, v);
    auto w2 = lila::ZerosLike(v);
    Apply(bonds, couplings, block, v, block, w2);
    REQUIRE(lila::close(w1, w2));

    auto evals_mat = lila::EigenvaluesSym(H);
    double e0_mat = evals_mat(0);
    double e0_app = E0Real(bonds, couplings, block);
    REQUIRE(lila::close(e0_mat, e0_app));
  }
}


TEST_CASE("spinhalf_apply", "[models][spinhalf]") {
  using namespace hydra::testcases::spinhalf;

  lila::Log.out("spinhalf_apply: Heisenberg chain apply test, J=1.0, N=2,..,6");
  for (int N = 2; N <= 6; ++N) {
    auto [bonds, couplings] = HBchain(N, 1.0);
    test_apply(bonds, couplings);
  }

  lila::Log.out("spinhalf_apply: Heisenberg alltoall apply test, N=2,..,6");
  for (int N = 2; N <= 6; ++N) {
    auto [bonds, couplings] = HB_alltoall(N);
    test_apply(bonds, couplings);
  }
}
