#include "../../catch.hpp"

#include <iostream>

#include "testcases_tj.h"
#include <hydra/all.h>
#include <lila/all.h>

using namespace hydra;

TEST_CASE("tJ_Apply", "[tj]") {
  using namespace hydra::testcases::tj;

  for (int N = 3; N <= 6; ++N) {
    lila::Log.out("TJModel: random all-to-all real apply=matrix, N={}", N);

    auto [bonds, cpls] = tj_alltoall(N);
    for (int nup = 0; nup <= N; ++nup)
      for (int ndn = 0; ndn <= N - nup; ++ndn) {

        auto block = tJ<uint32>(N, nup, ndn);
        auto H = MatrixReal(bonds, cpls, block, block);
        REQUIRE(lila::close(H, lila::Herm(H)));

        auto v = lila::Random<double>(block.size());
        auto w1 = lila::Mult(H, v);
        auto w2 = lila::ZerosLike(v);
        Apply(bonds, cpls, block, v, block, w2);
        REQUIRE(lila::close(w1, w2));

        auto evals_mat = lila::EigenvaluesSym(H);
        double e0_mat = evals_mat(0);
        double e0_app = E0Real(bonds, cpls, block);
        // lila::Log.out("nup: {}, ndn: {}, e0_mat: {}, e0_app: {}", nup, ndn,
        //              e0_mat, e0_app);
        REQUIRE(lila::close(e0_mat, e0_app));
      }
  }

  for (int N = 3; N <= 6; ++N) {
    lila::Log.out("TJModel: random all-to-all complex apply=matrix, N={}", N);

    auto [bonds, cpls] = tj_alltoall_complex(N);
    for (int nup = 0; nup <= N; ++nup)
      for (int ndn = 0; ndn <= N - nup; ++ndn) {
        auto block = tJ<uint32>(N, nup, ndn);
        auto H = MatrixCplx(bonds, cpls, block, block);
        REQUIRE(lila::close(H, lila::Herm(H)));

        auto v = lila::Random<complex>(block.size());
        auto w1 = lila::Mult(H, v);
        auto w2 = lila::ZerosLike(v);
        Apply(bonds, cpls, block, v, block, w2);
        REQUIRE(lila::close(w1, w2));

        auto evals_mat = lila::EigenvaluesSym(H);
        double e0_mat = evals_mat(0);
        double e0_app = E0Cplx(bonds, cpls, block);
        // lila::Log.out("nup: {}, ndn: {}, e0_mat: {}, e0_app: {}", nup, ndn,
        //              e0_mat, e0_app);
        REQUIRE(lila::close(e0_mat, e0_app));
      }
  }
}
