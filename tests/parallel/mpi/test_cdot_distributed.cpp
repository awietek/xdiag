#include <mpi.h>

#include <tests/catch.hpp>
#include <xdiag/parallel/mpi/allreduce.hpp>
#include <xdiag/parallel/mpi/cdot_distributed.hpp>

#include <xdiag/utils/print_macro.hpp>

using namespace xdiag;

template <class coeff_t>
void test_stable_dot(int size){
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  auto v = arma::Col<coeff_t>(size + rank, arma::fill::randu);
  auto w = arma::Col<coeff_t>(size + rank, arma::fill::randu);
  
  auto sdot = cdot_distributed(v, w);
  auto dot_proc = cdot(v, w);
  coeff_t dot;
  mpi::Allreduce(&dot_proc, &dot, 1, MPI_SUM, MPI_COMM_WORLD);
  REQUIRE(std::abs(dot - sdot) < 1e-12);    
}


TEST_CASE("cdot_distributed", "[mpi]") {
  Log("cdot_distributed test");
  for (int N = 2; N <= 6; ++N) {
    test_stable_dot<double>(N);
    test_stable_dot<complex>(N);
  }

}
