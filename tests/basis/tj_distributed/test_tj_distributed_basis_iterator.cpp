#include "../../catch.hpp"

#include <iostream>

#include "../../blocks/electron/testcases_electron.hpp"
#include <iostream>
#include <xdiag/blocks/tj_distributed.hpp>

using namespace xdiag;

void test_tj_distributed_basis_iterator(tJDistributed const &block) {
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  int64_t n_sites = block.n_sites();
  auto pprev = ProductState(n_sites);
  int64_t idx = 0;
  for (auto const &p : block) {
    // std::cout << mpi_rank << " " << to_string(p) << "\n";
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

TEST_CASE("tj_distributed_basis_iterator", "[basis]") {

  Log("tJDistributed Basis Iterator Np");
  for (int n_sites = 1; n_sites < 9; ++n_sites) {
    for (int n_up = 0; n_up <= n_sites; ++n_up) {
      for (int n_dn = 0; n_dn <= n_sites - n_up; ++n_dn) {
        auto block = tJDistributed(n_sites, n_up, n_dn);
        test_tj_distributed_basis_iterator(block);
      }
    }
  }
}
