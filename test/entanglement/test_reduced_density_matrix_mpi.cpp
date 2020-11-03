#include "../catch.hpp"

#include <hydra/all.h>
#include <iomanip>
#include <lila/allmpi.h>
#include <random>

#include "testcases_entanglement_entropy.h"

using namespace lila;
using namespace hydra;
using namespace hydra::entanglementtestcases;

template <class coeff_t>
void test_tJ_reduced_density_matrix_mpi(int N, BondList bondlist,
                                        Couplings couplings, qn_tj qn) {
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  // Check whether reduced density matrices agree with full ED
  for (int n_A_sites = 0; n_A_sites <= N; ++n_A_sites) {
    int n_B_sites = N - n_A_sites;

    for (number_t n_A_up = 0; n_A_up <= n_A_sites; ++n_A_up) {
      for (number_t n_A_dn = 0; n_A_dn <= n_A_sites; ++n_A_dn) {
        auto qn_A = qn_tj({n_A_up, n_A_dn});
        auto qn_B = qn - qn_A;

        if (valid(qn_A, n_A_sites) && valid(qn_B, n_B_sites)) {

          ///////////////////////////////////////////
          // Compute reduced density matrix from MPI
          auto H = TJModelMPI<coeff_t>(bondlist, couplings, qn);
          auto multiply_H = [&H](const VectorMPI<coeff_t> &v,
                                 VectorMPI<coeff_t> &w) {
            H.apply_hamiltonian(v, w);
          };
          // Create normal distributed random start state
          VectorMPI<coeff_t> startstate(H.local_dim());
          normal_dist_t<coeff_t> dist(0., 1.);
          normal_gen_t<coeff_t> gen(dist, 42 + mpi_rank);
          Random(startstate, gen, true);
          Normalize(startstate);

          // Run Lanczos
          auto res = LanczosEigenvectors(multiply_H, startstate, gen, false);
          auto gs = res.vectors[0];
          auto rho_A_mpi = ReducedDensityMatrix(H, n_A_sites, qn_A, gs);

          //////////////////////////////////////////////
          // Compute reduced density matrix from full ED
          auto Hfull = TJModel<coeff_t>(bondlist, couplings, qn);

          auto model = TJModel<coeff_t>(bondlist, couplings, qn);
          auto H2 = model.matrix();
          REQUIRE(lila::close(H2, lila::Herm(H2)));
          auto res2 = lila::EigenSym(H2);
          auto gs2 = res2.eigenvectors.col(0);
          auto rho_A_full = ReducedDensityMatrix(model, n_A_sites, qn_A, gs2);
          // if(mpi_rank==0)
          //   {
          //     LilaPrint(rho_A_full);
          //     LilaPrint(rho_A_mpi);
	  //     LilaPrint(rho_A_mpi -rho_A_full);
          //   }
          REQUIRE(lila::close(rho_A_full, rho_A_mpi, 1e-8, 1e-8));
        }
      }
    }
  }
}

TEST_CASE("ReducedDensityMatrixMPI", "[Entanglement]") {

  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  std::cout << std::setprecision(12);
  BondList bondlist;
  Couplings couplings;

  // Test entropy of t-J chains
  for (int N = 3; N <= 6; ++N) {
    if (mpi_rank == 0)
      std::cout << "@ ReducedDensityMatrixMPI tJ chain N=" << N << "\n";

    for (int n_up = 0; n_up <= N; ++n_up)
      for (int n_dn = 0; n_dn <= N; ++n_dn) {
        qn_tj qn = {n_up, n_dn};

        if (valid(qn, N)) {
          std::tie(bondlist, couplings) = tJchain_model(N);
          test_tJ_reduced_density_matrix_mpi<double>(N, bondlist, couplings,
                                                     qn);
        }
      }
  }

  // Test entropy of real random all-to-all tJ model
  for (int N = 3; N <= 6; ++N) {
    std::tie(bondlist, couplings) = tJrandom_model(N);
    if (mpi_rank == 0)
      std::cout << "@ ReducedDensityMatrixMPI tJ random (real) N=" << N << "\n";

    Couplings couplings_real;
    for (auto it : couplings)
      couplings_real[it.first] = lila::real(it.second);

    for (int n_up = 0; n_up <= N; ++n_up)
      for (int n_dn = 0; n_dn <= N; ++n_dn) {
        qn_tj qn = {n_up, n_dn};

        if (valid(qn, N)) {
          test_tJ_reduced_density_matrix_mpi<double>(N, bondlist,
                                                     couplings_real, qn);
        }
      }
  }

  // Test entropy of complex random all-to-all tJ model
  for (int N = 3; N <= 6; ++N) {
    std::tie(bondlist, couplings) = tJrandom_model(N);
    if (mpi_rank == 0)
      std::cout << "@ ReducedDensityMatrixMPI tJ random (cplx) N=" << N << "\n";

    for (int n_up = 0; n_up <= N; ++n_up)
      for (int n_dn = 0; n_dn <= N; ++n_dn) {
        qn_tj qn = {n_up, n_dn};

        if (valid(qn, N)) {
          double itensor_e0;
          std::vector<double> itensor_svns;
          std::tie(itensor_e0, itensor_svns) = tJrandom_e0_entropies(N, qn);
          test_tJ_reduced_density_matrix_mpi<complex>(N, bondlist, couplings,
                                                      qn);
        }
      }
  }
}
