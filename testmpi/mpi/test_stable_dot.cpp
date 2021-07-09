#include "../catch.hpp"

#include <hydra/allmpi.h>

using namespace hydra;

template <class coeff_t>
void test_stable_dot(int size){
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  auto v = lila::Random<coeff_t>(size + rank);
  auto w = lila::Random<coeff_t>(size + rank);
  
  auto sdot = mpi::StableDot(v, w);

  auto dot_proc = lila::Dot(v, w);
  coeff_t dot;
  mpi::Allreduce(&dot_proc, &dot, 1, MPI_SUM, MPI_COMM_WORLD);
  REQUIRE(lila::close(dot, sdot));    
}


TEST_CASE("stable_dot", "[mpi]") {

  LogMPI.out("StableDot test");

  for (int N = 2; N <= 6; ++N) {
    test_stable_dot<float>(N);
    test_stable_dot<double>(N);
    test_stable_dot<scomplex>(N);
    test_stable_dot<complex>(N);
  }

}
