// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"
#include <mpi.h>

#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/blocks/spinhalf_distributed.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/parallel/mpi/allreduce.hpp>

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
