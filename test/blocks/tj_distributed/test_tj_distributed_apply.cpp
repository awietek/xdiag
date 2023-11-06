#include "../../catch.hpp"

#include <iostream>

#include "../tj/testcases_tj.h"
#include <hydra/algorithms/sparse_diag.h>
#include <hydra/blocks/spinhalf/spinhalf_matrix.h>
#include <hydra/blocks/tj/tj_apply.h>
#include <hydra/blocks/tj/tj_matrix.h>
#include <hydra/blocks/tj_distributed/tj_distributed_apply.h>
#include <hydra/utils/close.h>
#include <hydra/utils/print_macro.h>

using namespace hydra;

TEST_CASE("tj_distributed_apply", "[tj_distributed]") try {
  using namespace hydra::testcases::tj;

  int N = 4;
  auto block = tJDistributed(N, N / 2, N / 2);
  BondList bonds;
  for (int i = 0; i < N; ++i) {
    bonds << Bond("ISING", "JZ", {i, (i + 1) % N});
  }
  bonds["JZ"] = 1.0;
  double e0 = eigval0(bonds, block);
  HydraPrint(e0);

  auto block2 = tJ(N, N / 2, N / 2);
  double e02 = eigval0(bonds, block2);
  HydraPrint(e02);
  HydraPrint(dim(block.basis()));
  HydraPrint(block2.dim());
  
} catch (std::exception const &e) {
  traceback(e);
}
