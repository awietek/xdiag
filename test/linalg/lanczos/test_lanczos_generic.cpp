#include "../../catch.hpp"

#include <hydra/all.h>

using namespace hydra;
using namespace lila;

template <class coeff_t>
void test_lanczos_generic(int n, int n_eigenvalue, real_t<coeff_t> precision,
                          int n_tests) {
  for (int random_seed : range<int>(n_tests)) {
    uniform_dist_t<coeff_t> fdist(-1., 1.);
    uniform_gen_t<coeff_t> fgen(fdist, random_seed);

    Matrix<coeff_t> A(n, n);
    Random(A, fgen);
    A += Herm(A);
    REQUIRE(close(A, Herm(A)));
    auto true_eigs = EigenvaluesSym(A);

    Vector<coeff_t> v0(n);
    Random(v0, fgen, false);

    // Define multipliers
    auto mult = [&A](const Vector<coeff_t> &v, Vector<coeff_t> &w) {
      w = Mult(A, v);
    };
    auto dot = [](lila::Vector<coeff_t> const &v,
                  lila::Vector<coeff_t> const &w) -> coeff_t {
      return lila::Dot(v, w);
    };
    auto converged = [n_eigenvalue, precision](Tmatrix const &tmat) -> bool {
      return ConvergedEigenvalues(tmat, n_eigenvalue, precision);
    };
    // Lanczos eigenvalue computation with precision on eigenvalues
    auto [tmat, _] = LanczosGeneric(mult, v0, dot, converged);
    auto [teigs, tevecs] = tmat.eigen();
    REQUIRE(close(teigs(0), true_eigs(0)));

    // Rerun for eigenvectors, check whether energy agrees
    Random(v0, fgen, false);
    auto [tmat2, evec] =
        LanczosGeneric(mult, v0, dot, converged, tevecs.col(n_eigenvalue));
    auto &v = evec;
    auto Hv = v;
    mult(v, Hv);
    auto e = dot(v, Hv);

    REQUIRE(close(real(e), true_eigs(n_eigenvalue)));
    REQUIRE(close(imag(e), 0.0));
    REQUIRE(close(Norm(v), 1.0));

  }
}

TEST_CASE("lanczos_generic", "[lanczos]") {
  test_lanczos_generic<double>(100, 0, 1e-12, 2);
  test_lanczos_generic<complex>(100, 0, 1e-12, 2);
  test_lanczos_generic<double>(200, 1, 1e-12, 2);
  test_lanczos_generic<complex>(200, 1, 1e-12, 2);
  test_lanczos_generic<double>(400, 2, 1e-12, 2);
  test_lanczos_generic<complex>(400, 2, 1e-12, 2);
}
