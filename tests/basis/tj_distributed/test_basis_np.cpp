// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include <mpi.h>

#include <xdiag/basis/tj_distributed/basis_np.hpp>
#include <xdiag/bits/bitops.hpp>
#include <xdiag/combinatorics/binomial.hpp>
#include <xdiag/parallel/mpi/allreduce.hpp>
#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/utils/logger.hpp>

template <typename bit_t, typename coeff_t>
void test_tj_distributed_basis_np_transpose() {
  using namespace xdiag;

  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  for (int nsites = 1; nsites <= 10; ++nsites) {
    for (int nup = 0; nup <= nsites; ++nup) {
      for (int ndn = 0; ndn < nsites - nup; ++ndn) {

        auto basis = basis::tj_distributed::BasisNp<bit_t>(nsites, nup, ndn);
	int64_t buffer_size = std::max(basis.size(), basis.size_transpose());
	mpi::buffer.reserve<coeff_t>(buffer_size);

        arma::Col<coeff_t> v(basis.size(), arma::fill::randu);
        arma::Col<coeff_t> w(basis.size_transpose(), arma::fill::randu);
        arma::Col<coeff_t> w2(basis.size_transpose(), arma::fill::randu);

        basis.transpose(v.memptr());
        std::copy(mpi::buffer.send<coeff_t>(),
                  mpi::buffer.send<coeff_t>() + basis.size_transpose(),
                  w.memptr());

        basis.transpose(v.memptr(), w2.memptr());
        REQUIRE(isapprox(w, w2));

        ////// DEBUG PRINT
        // for (int rank = 0; rank < mpi_size; ++rank) {
        //   if (rank == mpi_rank) {
        //     std::cout << "[" << rank << "]:\n";
        //     int64_t idx_up = 0;
        //     int64_t idx = 0;
        //     for (bit_t up : basis.my_ups()) {
        //       for (bit_t dn : basis.my_dns_for_ups(idx_up)) {
        //         std::cout << up << " " << dn << "  v(" << idx << ")=" <<
        //         v(idx)
        //                   << "\n";
        //         ++idx;
        //       }
        //       ++idx_up;
        //     }
        //   }
        //   MPI_Barrier(MPI_COMM_WORLD);
        // }
        // if (mpi_rank == 0) {
        //   std::cout << "\n";
        // }

        // for (int rank = 0; rank < mpi_size; ++rank) {
        //   if (rank == mpi_rank) {
        //     std::cout << "[" << rank << "]:\n";
        //     int64_t idx_dn = 0;
        //     int64_t idx = 0;
        //     for (bit_t dn : basis.my_dns()) {
        //       for (bit_t up : basis.my_ups_for_dns(idx_dn)) {
        //         std::cout << up << " " << dn << "  w(" << idx << ")=" <<
        //         w(idx)
        //                   << "\n";
        //         ++idx;
        //       }
        //       ++idx_dn;
        //     }
        //   }
        //   MPI_Barrier(MPI_COMM_WORLD);
        // }
        // if (mpi_rank == 0) {
        //   std::cout << "------------------------------\n";
        //   std::cout << "------------------------------\n";
        //   std::cout << "------------------------------\n";
        // }

        arma::Col<coeff_t> v2(basis.size(), arma::fill::randu);
        arma::Col<coeff_t> v3(basis.size(), arma::fill::randu);

        basis.transpose_r(w.memptr());
        std::copy(mpi::buffer.send<coeff_t>(),
                  mpi::buffer.send<coeff_t>() + basis.size(), v2.memptr());

        basis.transpose_r(w.memptr(), v3.memptr());
        REQUIRE(isapprox(v2, v3));

        ////////// DEBUG PRINT
        // for (int rank = 0; rank < mpi_size; ++rank) {
        //   if (rank == mpi_rank) {
        //     std::cout << "[" << rank << "]:\n";
        //     int64_t idx_dn = 0;
        //     int64_t idx = 0;
        //     for (bit_t dn : basis.my_dns()) {
        //       for (bit_t up : basis.my_ups_for_dns(idx_dn)) {
        //         std::cout << up << " " << dn << "  w(" << idx << ")=" <<
        //         w(idx)
        //                   << "\n";
        //         ++idx;
        //       }
        //       ++idx_dn;
        //     }
        //   }
        //   MPI_Barrier(MPI_COMM_WORLD);
        // }
        // if (mpi_rank == 0) {
        //   std::cout << "\n";
        // }

        // for (int rank = 0; rank < mpi_size; ++rank) {
        //   if (rank == mpi_rank) {
        //     std::cout << "[" << rank << "]:\n";
        //     int64_t idx_up = 0;
        //     int64_t idx = 0;
        //     for (bit_t up : basis.my_ups()) {
        //       for (bit_t dn : basis.my_dns_for_ups(idx_up)) {
        //         std::cout << up << " " << dn << "  v(" << idx << ")=" <<
        //         v2(idx)
        //                   << "\n";
        //         ++idx;
        //       }
        //       ++idx_up;
        //     }
        //   }
        //   MPI_Barrier(MPI_COMM_WORLD);
        // }

        REQUIRE(isapprox(v2, v));
      }
    }
  }
}

template <typename bit_t> void test_tj_distributed_basis_np() {
  using namespace xdiag;

  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  arma::arma_rng::set_seed(mpi_rank);

  for (int nsites = 0; nsites <= 6; ++nsites) {
    for (int nup = 0; nup <= nsites; ++nup) {
      for (int ndn = 0; ndn <= nsites - nup; ++ndn) {
        // int nsites=16;
        // int nup = 7;
        // int ndn = 6;
        auto basis = basis::tj_distributed::BasisNp<bit_t>(nsites, nup, ndn);

        // ups/ dns order
        bit_t ups_before = 0;
        int64_t size_local = 0;
        for (int64_t idx_ups = 0; idx_ups < (int64_t)basis.my_ups().size();
             ++idx_ups) {
          bit_t ups = basis.my_ups()[idx_ups];
          REQUIRE(xdiag::bits::popcnt(ups) == nup);
          if (idx_ups != 0) {
            REQUIRE(ups > ups_before);
          }
          auto dnss = basis.my_dns_for_ups(idx_ups);

          int64_t idx_dns = 0;
          bit_t dns_before = 0;
          for (auto dns : dnss) {
            REQUIRE(xdiag::bits::popcnt(dns) == ndn);
            if (idx_dns != 0) {
              REQUIRE(dns > dns_before);
            }
            ++idx_dns;
            ++size_local;
          }
          REQUIRE(idx_dns == combinatorics::binomial(nsites - nup, ndn));
          ups_before = ups;
        }
        REQUIRE(size_local == basis.size());
        int64_t total_size = 0;
        mpi::Allreduce(&size_local, &total_size, 1, MPI_SUM, MPI_COMM_WORLD);
        REQUIRE(total_size == basis.dim());

        // dns/ ups order
        bit_t dns_before = 0;
        int64_t size_local_transpose = 0;
        for (int64_t idx_dns = 0; idx_dns < (int64_t)basis.my_dns().size();
             ++idx_dns) {
          bit_t dns = basis.my_dns()[idx_dns];
          REQUIRE(xdiag::bits::popcnt(dns) == ndn);
          if (idx_dns != 0) {
            REQUIRE(dns > dns_before);
          }
          auto upss = basis.my_ups_for_dns(idx_dns);

          int64_t idx_ups = 0;
          bit_t ups_before = 0;
          for (auto ups : upss) {
            REQUIRE(xdiag::bits::popcnt(ups) == nup);
            if (idx_ups != 0) {
              REQUIRE(ups > ups_before);
            }
            ++idx_ups;
            ++size_local_transpose;
          }
          REQUIRE(idx_ups == combinatorics::binomial(nsites - ndn, nup));
          dns_before = dns;
        }
        REQUIRE(size_local_transpose == basis.size_transpose());
        total_size = 0;
        mpi::Allreduce(&size_local_transpose, &total_size, 1, MPI_SUM,
                       MPI_COMM_WORLD);
        REQUIRE(total_size == basis.dim());

        // Log("{} {} {} {} {}", nsites, nup, ndn, basis.size_max(),
        //     basis.size_min());
      }
    }
  }
}

TEST_CASE("tj_distributed_basis_np", "[tj_distributed]") {
  using namespace xdiag;

  Log("tj_distributed_basis_np transpose test (uint32_t, double)");
  test_tj_distributed_basis_np_transpose<uint32_t, double>();
  Log("tj_distributed_basis_np transpose test (uint64_t, double)");
  test_tj_distributed_basis_np_transpose<uint64_t, double>();

  Log("tj_distributed_basis_np transpose test (uint32_t, complex)");
  test_tj_distributed_basis_np_transpose<uint32_t, complex>();
  Log("tj_distributed_basis_np transpose test (uint64_t, complex)");
  test_tj_distributed_basis_np_transpose<uint64_t, complex>();

  Log("tj_distributed_basis_np test (uint32_t)");
  test_tj_distributed_basis_np<uint32_t>();
  Log("tj_distributed_basis_np test (uint64_t)");
  test_tj_distributed_basis_np<uint64_t>();
}
