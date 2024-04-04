#include "../catch.hpp"

#include <xdiag/all.h>
#include <iomanip>
#include <lila/all.h>
#include <random>

#include "testcases_entanglement_entropy.h"

using namespace lila;
using namespace xdiag;
using namespace xdiag::entanglementtestcases;

template <class coeff_t>
void test_tJ_entanglement_entropy(int N, BondList bondlist, Couplings couplings,
                                  double itensor_e0,
                                  std::vector<double> itensor_svns, qn_tj qn) {
  // std::cout << String(qn) << "\n";

  auto model = TJModel<coeff_t>(bondlist, couplings, qn);
  auto H = model.matrix();

  REQUIRE(lila::close(H, lila::Herm(H)));
  
  auto res = lila::EigenSym(H);
  auto e0 = res.eigenvalues(0);
  auto gs = res.eigenvectors.col(0);

  CHECK(std::abs(e0 - itensor_e0) < 1e-10);

  // if (std::abs(e0 - itensor_e0) > 1e-10)
  //   LilaPrint(res.eigenvalues);

  // std::cout << std::setprecision(16);
  // std::cout << "  ED      E0 " << e0 << "\n";
  // std::cout << "  ITensor E0 " << itensor_e0 << "\n";

  for (int b = 1; b < N; ++b) {
    auto svn = EntanglementEntropy(model, gs, b);
    // std::cout << "  ED      b=" << b << ", SvN = " << svn << "\n";
    // std::cout << "  ITensor b=" << b << ", SvN = " << itensor_svns[b - 1]
    //           << "\n";
    CHECK(std::abs(svn - itensor_svns[b - 1]) < 1e-6);
  }
  // std::cout << "\n";
}

TEST_CASE("EntanglementEntropy", "[Entanglement]") {

  std::cout << std::setprecision(12);
  BondList bondlist;
  Couplings couplings;

  // Test entropy of t-J chains
  for (int N = 3; N <= 6; ++N) {
    std::cout << "@ EntanglementEntropy tJ chain N=" << N << "\n";

    for (int n_up = 0; n_up <= N; ++n_up)
      for (int n_dn = 0; n_dn <= N; ++n_dn) {
        qn_tj qn = {n_up, n_dn};

        if (valid(qn, N)) {
          double itensor_e0;
          std::vector<double> itensor_svns;
          std::tie(itensor_e0, itensor_svns) = tJchain_e0_entropies(N, qn);
          std::tie(bondlist, couplings) = tJchain_model(N);
          test_tJ_entanglement_entropy<double>(N, bondlist, couplings,
                                               itensor_e0, itensor_svns, qn);
        }
      }
  }

  // Test entropy of real random all-to-all tJ model
  for (int N = 3; N <= 6; ++N) {
    std::tie(bondlist, couplings) = tJrandom_model(N);
    std::cout << "@ EntanglementEntropy tJ random (real) N=" << N << "\n";

    Couplings couplings_real;
    for (auto it : couplings)
      couplings_real[it.first] = lila::real(it.second);

    for (int n_up = 0; n_up <= N; ++n_up)
      for (int n_dn = 0; n_dn <= N; ++n_dn) {
        qn_tj qn = {n_up, n_dn};

        if (valid(qn, N)) {
          double itensor_e0;
          std::vector<double> itensor_svns;
          std::tie(itensor_e0, itensor_svns) =
              tJrandom_e0_entropies(N, qn, true);
          test_tJ_entanglement_entropy<double>(N, bondlist, couplings_real,
                                               itensor_e0, itensor_svns, qn);
        }
      }
  }

  // Test entropy of complex random all-to-all tJ model
  for (int N = 3; N <= 6; ++N) {
    std::tie(bondlist, couplings) = tJrandom_model(N);
    std::cout << "@ EntanglementEntropy tJ random (cplx) N=" << N << "\n";

    for (int n_up = 0; n_up <= N; ++n_up)
      for (int n_dn = 0; n_dn <= N; ++n_dn) {
        qn_tj qn = {n_up, n_dn};

        if (valid(qn, N)) {
          double itensor_e0;
          std::vector<double> itensor_svns;
          std::tie(itensor_e0, itensor_svns) = tJrandom_e0_entropies(N, qn);
          test_tJ_entanglement_entropy<complex>(N, bondlist, couplings,
                                                itensor_e0, itensor_svns,
                                                qn);
        }
      }
  }
}
