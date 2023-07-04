#include "../../catch.hpp"

#include <iostream>

#include <hydra/algorithms/lanczos/lanczos_pro.h>
#include <hydra/utils/print_macro.h>

TEST_CASE("lanczos_pro", "[lanczos]") {
  using namespace hydra;
  using namespace arma;

  Log("testing lanczos_pro");
  
  int N = 100;
  mat B = randn(N, N);
  mat A = B + B.t();

  auto mult = [&A](vec const &v, vec &w) { w = A * v; };

  vec v0 = randn(N);
  auto conv = []() { return true; };

  auto res = lanczos_pro(mult, v0, conv);

  mat V = res.V;
  HydraPrint(mat(V.t() * V));
  HydraPrint(res.omegas);  
}
