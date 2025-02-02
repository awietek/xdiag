#include "../../catch.hpp"

#include <iostream>

#include "../../blocks/electron/testcases_electron.hpp"
#include <iostream>
#include <xdiag/blocks/spinhalf_distributed.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;

void test_spinhalf_distributed_basis_iterator(
    SpinhalfDistributed const &block) {
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  int64_t nsites = block.nsites();
  auto pprev = ProductState(nsites);
  int64_t idx = 0;
  for (auto const &p : block) {
    REQUIRE(p != pprev);
    int64_t idx2 = block.index(p);
    // std::cout << mpi_rank << " " << to_string(p) << " " << idx << " " << idx2
    //           << "\n";
    REQUIRE(idx == idx2);
    pprev = p;
    ++idx;
  }
  REQUIRE(idx == block.size());
}

TEST_CASE("spinhalf_distributed_basis_iterator", "[basis]") {

  Log("SpinhalfDistributed Basis Iterator Sz");
  for (int nsites = 1; nsites < 9; ++nsites) {
    for (int nup = 0; nup <= nsites; ++nup) {
      auto block = SpinhalfDistributed(nsites, nup);
      test_spinhalf_distributed_basis_iterator(block);
    }
  }
}
