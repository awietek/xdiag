#include "../../catch.hpp"

#include <iostream>

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/apply.hpp>
#include <xdiag/algorithms/lanczos/lanczos_pro.hpp>
#include <xdiag/algorithms/sparse_diag.hpp>
#include <xdiag/io/read.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/create_state.hpp>
#include <xdiag/states/random_state.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/representation.hpp>
#include <xdiag/algebra/isapprox.hpp>

TEST_CASE("lanczos_pro", "[lanczos]") {
  using namespace xdiag;
  using namespace arma;

  Log.set_verbosity(0);
  Log("testing lanczos_pro");

  {
    Log("  random Hermitian orthogonal");
    int N = 2000;
    mat B = randn(N, N);
    mat A = B + B.t();

    auto mult = [&A](vec const &v, vec &w) { w = A * v; };

    vec v0 = randn(N);
    auto conv = [](Tmatrix) { return false; };

    double ortho_level = 1e-8;

    int niter = 500;
    auto res = lanczos_pro(mult, v0, conv, niter, ortho_level, 1e-7, true);
    // XDIAG_SHOW(res.num_reorthogonalizations);
    mat V = res.V;
    double total_orthogonality = norm(V.t() * V - eye(niter, niter));
    REQUIRE(total_orthogonality < niter * ortho_level);
  }

  {
    Log("  random Hermitian exact deflation, exact V A Vt decomp");
    int N = 567;
    mat B = randn(N, N);
    mat A = B + B.t();

    auto mult = [&A](vec const &v, vec &w) { w = A * v; };

    vec v0 = randn(N);
    auto conv = [](Tmatrix) { return false; };

    double ortho_level = 1e-8;

    int niter = 1111;
    auto res = lanczos_pro(mult, v0, conv, niter, ortho_level, 1e-7, true);
    // XDIAG_SHOW(res.num_reorthogonalizations);

    REQUIRE(res.num_iterations == N); // exact deflation

    mat V = res.V;
    double total_orthogonality = norm(V.t() * V - eye(N, N));
    REQUIRE(total_orthogonality < N * ortho_level); // correct ortho levels

    vec eigs;
    eig_sym(eigs, A);
    vec eigs_lcs = res.eigenvalues;
    // XDIAG_SHOW(eigs);
    // XDIAG_SHOW(eigs_lcs);
    REQUIRE(isapprox(eigs, eigs_lcs)); // exact eigenvalue solution

    mat tmat = res.tmat;
    mat A2 = mat(V * tmat * V.t());
    REQUIRE(norm(A - A2) < 2 * N * ortho_level);
  }

  {
    Log("  random Hermitian orthogonal (complex)");
    int N = 2000;
    cx_mat B(N, N, fill::randn);
    cx_mat A = B + B.t();

    auto mult = [&A](cx_vec const &v, cx_vec &w) { w = A * v; };

    cx_vec v0(N, fill::randn);
    auto conv = [](Tmatrix) { return false; };

    double ortho_level = 1e-8;

    int niter = 500;
    auto res = lanczos_pro(mult, v0, conv, niter, ortho_level, 1e-7, true);
    // XDIAG_SHOW(res.num_reorthogonalizations);

    cx_mat V = res.V;
    double total_orthogonality = norm(V.t() * V - eye(niter, niter));
    REQUIRE(total_orthogonality < niter * ortho_level);
  }

  {
    Log("  random Hermitian exact deflation, exact V A Vt decomp (complex)");
    int N = 567;
    cx_mat B(N, N, fill::randn);
    cx_mat A = B + B.t();

    auto mult = [&A](cx_vec const &v, cx_vec &w) { w = A * v; };

    cx_vec v0(N, fill::randn);
    auto conv = [](Tmatrix) { return false; };

    double ortho_level = 1e-8;

    int niter = 1111;
    auto res = lanczos_pro(mult, v0, conv, niter, ortho_level, 1e-7, true);
    // XDIAG_SHOW(res.num_reorthogonalizations);

    REQUIRE(res.num_iterations == N); // exact deflation

    cx_mat V = res.V;
    double total_orthogonality = norm(V.t() * V - eye(N, N));
    REQUIRE(total_orthogonality < 2 * N * ortho_level); // correct ortho levels

    vec eigs;
    eig_sym(eigs, A);
    vec eigs_lcs = res.eigenvalues;
    // XDIAG_SHOW(eigs);
    // XDIAG_SHOW(eigs_lcs);
    REQUIRE(isapprox(eigs, eigs_lcs)); // exact eigenvalue solution

    cx_mat tmatc(res.tmat, mat(N, N, fill::zeros));
    cx_mat A2 = cx_mat(V * tmatc * V.t());
    REQUIRE(norm(A - A2) < 2 * N * ortho_level);
  }

  {
    Log("  shastry sutherland N=16 ortho levels");

    std::string lfilename =
        XDIAG_DIRECTORY "/misc/data/shastry.16.HB.J.Jd.fsl.toml";
    auto lfile = FileToml(lfilename);
    auto ops = lfile["Interactions"].as<OpSum>();
    int nsites = 16;
    ops["J"] = 0.63;
    ops["Jd"] = 1.00;

    std::vector<std::string> irrep_names = {"Gamma.C1.A", "M.C1.A", "X0.C1.A",
                                            "X1.C1.A"};

    for (int nup = 0; nup <= 6; ++nup) {
      for (std::string k : irrep_names) {
        auto irrep = read_representation(lfile, k);
        auto block = Spinhalf(nsites, nup, irrep);

        Log("   nup: {} k: {}", nup, k);
        // XDIAG_SHOW(block);

        auto rstate = random_state(block, false);
        auto mult = [&ops, &block](cx_vec const &v, cx_vec &w) {
          apply(ops, block, v, block, w);
        };
        auto conv = [](Tmatrix) { return false; };
        int n_iter = 500;

        if (block.size() < n_iter) {
          continue;
        }

        cx_vec v0 = rstate.vectorC();

        double ortho_level = (n_iter > block.size()) ? 0 : 1e-8;
        auto res =
            lanczos_pro(mult, v0, conv, n_iter, ortho_level, 1e-7, false);

        // XDIAG_SHOW(res.num_reorthogonalizations);
        int n_iter2 = res.num_iterations;

        cx_mat V = res.V;
        // XDIAG_SHOW(res.num_iterations);
        // XDIAG_SHOW(cx_mat(V.t() * V ));
        double total_orthogonality = norm(V.t() * V - eye(n_iter2, n_iter2));
        if (n_iter > block.size()) {
          CHECK(total_orthogonality < 1e-6);
        } else {
          CHECK(total_orthogonality < 2 * n_iter * ortho_level);
        }
      }
    }
  }

  {
    Log("  shastry sutherland N=20 ortho levels");

    std::string lfilename =
        XDIAG_DIRECTORY "/misc/data/shastry.20.HB.J.Jd.fsl.toml";
    auto lfile = FileToml(lfilename);
    auto ops = lfile["Interactions"].as<OpSum>();
    int nsites = 20;
    ops["J"] = 0.63;
    ops["Jd"] = 1.00;

    std::vector<std::string> irrep_names = {"Gamma.C1.A"};

    for (int nup = 0; nup <= 6; ++nup) {
      for (std::string k : irrep_names) {
        auto irrep = read_representation(lfile, k);
        auto block = Spinhalf(nsites, nup, irrep);

        Log("   nup: {} k: {}", nup, k);
        // XDIAG_SHOW(block);

        auto rstate = random_state(block, false);
        auto mult = [&ops, &block](cx_vec const &v, cx_vec &w) {
          apply(ops, block, v, block, w);
        };
        auto conv = [](Tmatrix) { return false; };
        int n_iter = 500;

        if (block.size() < n_iter) {
          continue;
        }

        cx_vec v0 = rstate.vectorC();

        double ortho_level = (n_iter > block.size()) ? 0 : 1e-8;
        auto res =
            lanczos_pro(mult, v0, conv, n_iter, ortho_level, 1e-7, false);

        // XDIAG_SHOW(res.num_reorthogonalizations);
        int n_iter2 = res.num_iterations;
        cx_mat V = res.V;
        double total_orthogonality = norm(V.t() * V - eye(n_iter2, n_iter2));
        if (n_iter > block.size()) {
          CHECK(total_orthogonality < 1e-6);
        } else {
          CHECK(total_orthogonality < n_iter * ortho_level);
        }
      }
    }
  }
}
