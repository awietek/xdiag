#include "../../catch.hpp"
#include <mpi.h>

#include <hydra/basis/tj_distributed/basis_np.h>
#include <hydra/bits/bitops.h>
#include <hydra/combinatorics/binomial.h>
#include <hydra/parallel/mpi/allreduce.h>
#include <hydra/utils/print_macro.h>

TEST_CASE("tj_distributed_basis_np", "[tj_distributed]") {
  using namespace hydra;
  using hydra::bits::popcnt;
  using bit_t = uint32_t;

  Log("tj_distributed_basis_np test");
  for (int n_sites = 0; n_sites <= 6; ++n_sites) {
    for (int n_up = 0; n_up <= n_sites; ++n_up) {
      for (int n_dn = 0; n_dn <= n_sites - n_up; ++n_dn) {
        // int n_sites=16;
        // int n_up = 7;
        // int n_dn = 6;
        auto basis = basis::tj_distributed::BasisNp<bit_t>(n_sites, n_up, n_dn);

        // ups/ dns order
        bit_t ups_before = 0;
        int64_t size_local = 0;
        for (int64_t idx_ups = 0; idx_ups < (int64_t)basis.my_ups().size();
             ++idx_ups) {
          bit_t ups = basis.my_ups()[idx_ups];
          REQUIRE(popcnt(ups) == n_up);
          if (idx_ups != 0) {
            REQUIRE(ups > ups_before);
          }
          auto dnss = basis.my_dns_for_ups(idx_ups);

          int64_t idx_dns = 0;
          bit_t dns_before = 0;
          for (auto dns : dnss) {
            REQUIRE(popcnt(dns) == n_dn);
            if (idx_dns != 0) {
              REQUIRE(dns > dns_before);
            }
            ++idx_dns;
            ++size_local;
          }
          REQUIRE(idx_dns == combinatorics::binomial(n_sites - n_up, n_dn));
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
          REQUIRE(popcnt(dns) == n_dn);
          if (idx_dns != 0) {
            REQUIRE(dns > dns_before);
          }
          auto upss = basis.my_ups_for_dns(idx_dns);

          int64_t idx_ups = 0;
          bit_t ups_before = 0;
          for (auto ups : upss) {
            REQUIRE(popcnt(ups) == n_up);
            if (idx_ups != 0) {
              REQUIRE(ups > ups_before);
            }
            ++idx_ups;
            ++size_local_transpose;
          }
          REQUIRE(idx_ups == combinatorics::binomial(n_sites - n_dn, n_up));
          dns_before = dns;
        }
        REQUIRE(size_local_transpose == basis.size_transpose());
        total_size = 0;
        mpi::Allreduce(&size_local_transpose, &total_size, 1, MPI_SUM, MPI_COMM_WORLD);
        REQUIRE(total_size == basis.dim());

        // Log("{} {} {} {} {}", n_sites, n_up, n_dn, basis.size_max(),
        //     basis.size_min());
      }
    }
  }
}
