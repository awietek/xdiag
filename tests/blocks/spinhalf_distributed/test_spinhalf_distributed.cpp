#include "../../catch.hpp"
#include <mpi.h>

#include <xdiag/all.hpp>

using namespace xdiag;

void test_spinhalf_distributed(int nsites) {
  for (int nup = 0; nup <= nsites; ++nup) {
    auto block = Spinhalf(nsites, nup);
    auto block_mpi = SpinhalfDistributed(nsites, nup);

    auto mysiz = block_mpi.size();
    int64_t mysize = 0;

    mpi::Allreduce(&mysiz, &mysize, 1, MPI_SUM, MPI_COMM_WORLD);

    REQUIRE(block_mpi.dim() == block.size());
    REQUIRE(block_mpi.dim() == mysize);
    REQUIRE(block_mpi.dim() == combinatorics::binomial(nsites, nup));
  }
}

TEST_CASE("spinhalf_distributed", "[spinhalf_distributed]") {

  Log("SpinhalfDistributed test");
  for (int N = 0; N <= 6; ++N) {
    test_spinhalf_distributed(N);
  }
}
