#include "../../catch.hpp"

#include <hydra/allmpi.h>

using namespace hydra;

template <class bit_t>
void test_spinhalf_mpi(int n_sites){
  for (int nup = 0; nup <= n_sites; ++nup) {
    auto block = Spinhalf<uint32>(n_sites, nup);
    auto block_mpi = SpinhalfMPI<uint32>(n_sites, nup);

    auto mysiz = block_mpi.size();
    lila::size_type mysize = 0;

    mpi::Allreduce(&mysiz, &mysize, 1, MPI_SUM, MPI_COMM_WORLD);

    REQUIRE(block_mpi.dim() == block.dim());
    REQUIRE(block_mpi.dim() == mysize);
    REQUIRE(block_mpi.dim() == combinatorics::binomial(n_sites, nup));
  }
}


TEST_CASE("spinhalf_mpi", "[spinhalf_mpi]") {

  LogMPI.out("SpinhalfMPI test");

  for (int N = 2; N <= 6; ++N) {
    test_spinhalf_mpi<uint32>(N);
  }

}
