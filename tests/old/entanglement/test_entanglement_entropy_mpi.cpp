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
void test_tJ_entanglement_entropy_mpi(int N, BondList bondlist,
                                      Couplings couplings, double itensor_e0,
                                      std::vector<double> itensor_svns,
                                      qn_tj qn) {
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  
  auto model = TJModelMPI<coeff_t>(bondlist, couplings, qn);
  auto multiply_H = [&model](const VectorMPI<coeff_t> &v,
                             VectorMPI<coeff_t> &w) {
    model.apply_hamiltonian(v, w);
  };
  // Create normal distributed random start state
  VectorMPI<coeff_t> startstate(model.local_dim());
  normal_dist_t<coeff_t> dist(0., 1.);
  normal_gen_t<coeff_t> gen(dist, 42 + mpi_rank);
  Random(startstate, gen, true);
  Normalize(startstate);

  // Run Lanczos
  auto res = LanczosEigenvectors(multiply_H, startstate, gen, false);
  auto gs = res.vectors[0];
  auto e0 = res.eigenvalues(0);

  CHECK(std::abs(e0 - itensor_e0) < 1e-10);

  for (int b = 1; b < N; ++b) {
    auto svn = EntanglementEntropy(model, gs, b);
    CHECK(std::abs(svn - itensor_svns[b - 1]) < 1e-6);
  }
}

TEST_CASE("EntanglementEntropyMPI", "[Entanglement]") {

  std::cout << std::setprecision(12);
  BondList bondlist;
  Couplings couplings;

  // Test entropy of t-J chains
  for (int N = 3; N <= 6; ++N) {
    std::cout << "@ EntanglementEntropyMPI tJ chain N="
              << N << "\n";

    for (int nup = 0; nup <= N; ++nup)
      for (int ndn = 0; ndn <= N; ++ndn) {
        qn_tj qn = {nup, ndn};

        if (valid(qn, N)) {
          double itensor_e0;
          std::vector<double> itensor_svns;
          std::tie(itensor_e0, itensor_svns) = tJchain_e0_entropies(N, qn);
          std::tie(bondlist, couplings) = tJchain_model(N);
          test_tJ_entanglement_entropy_mpi<double>(
              N, bondlist, couplings, itensor_e0, itensor_svns, qn);
        }
      }
  }

  // Test entropy of real random all-to-all tJ model
  for (int N = 3; N <= 6; ++N) {
    std::tie(bondlist, couplings) = tJrandom_model(N);
    std::cout
        << "@ EntanglementEntropyMPI tJ random (real) N="
        << N << "\n";

    Couplings couplings_real;
    for (auto it : couplings)
      couplings_real[it.first] = lila::real(it.second);

    for (int nup = 0; nup <= N; ++nup)
      for (int ndn = 0; ndn <= N; ++ndn) {
        qn_tj qn = {nup, ndn};

        if (valid(qn, N)) {
          double itensor_e0;
          std::vector<double> itensor_svns;
          std::tie(itensor_e0, itensor_svns) =
              tJrandom_e0_entropies(N, qn, true);
          test_tJ_entanglement_entropy_mpi<double>(
              N, bondlist, couplings_real, itensor_e0, itensor_svns, qn);
        }
      }
  }

  // Test entropy of complex random all-to-all tJ model
  for (int N = 3; N <= 6; ++N) {
    std::tie(bondlist, couplings) = tJrandom_model(N);
    std::cout
        << "@ EntanglementEntropyMPI tJ random (cplx) N="
        << N << "\n";

    for (int nup = 0; nup <= N; ++nup)
      for (int ndn = 0; ndn <= N; ++ndn) {
        qn_tj qn = {nup, ndn};

        if (valid(qn, N)) {
          double itensor_e0;
          std::vector<double> itensor_svns;
          std::tie(itensor_e0, itensor_svns) = tJrandom_e0_entropies(N, qn);
          test_tJ_entanglement_entropy_mpi<complex>(
              N, bondlist, couplings, itensor_e0, itensor_svns, qn);
        }
      }
  }
}
