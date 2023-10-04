#include <mpi.h>
#include "../../catch.hpp"

#include <hydra/basis/tj_distributed/basis_np.h>
#include <hydra/utils/print_macro.h>


TEST_CASE("tj_distributed_basis_np", "[tj_distributed]") {
  using namespace hydra;

  Log("tj_distributed_basis_np test");
  auto basis = basis::tj_distributed::BasisNp<uint32_t>(6, 3, 2);
  HydraPrint(basis.size());
}
