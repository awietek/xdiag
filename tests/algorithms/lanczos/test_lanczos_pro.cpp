#include "../../catch.hpp"

#include <iostream>

#include <hydra/algebra/algebra.h>
#include <hydra/algebra/apply.h>

#include <hydra/algorithms/lanczos/lanczos_pro.h>
#include <hydra/algorithms/sparse_diag.h>
#include <hydra/io/file_toml.h>
#include <hydra/operators/bondlist.h>
#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/representation.h>
#include <hydra/utils/close.h>
#include <hydra/utils/print_macro.h>

TEST_CASE("lanczos_pro", "[lanczos]") {
  using namespace hydra;
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
    // HydraPrint(res.num_reorthogonalizations);
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
    // HydraPrint(res.num_reorthogonalizations);

    REQUIRE(res.num_iterations == N); // exact deflation

    mat V = res.V;
    double total_orthogonality = norm(V.t() * V - eye(N, N));
    REQUIRE(total_orthogonality < N * ortho_level); // correct ortho levels

    vec eigs;
    eig_sym(eigs, A);
    vec eigs_lcs = res.eigenvalues;
    // HydraPrint(eigs);
    // HydraPrint(eigs_lcs);
    REQUIRE(close(eigs, eigs_lcs)); // exact eigenvalue solution

    mat tmat = res.tmat;
    mat A2 = mat(V * tmat * V.t());
    REQUIRE(norm(A - A2) < N * ortho_level);
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
    // HydraPrint(res.num_reorthogonalizations);

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
    // HydraPrint(res.num_reorthogonalizations);

    REQUIRE(res.num_iterations == N); // exact deflation

    cx_mat V = res.V;
    double total_orthogonality = norm(V.t() * V - eye(N, N));
    REQUIRE(total_orthogonality < N * ortho_level); // correct ortho levels

    vec eigs;
    eig_sym(eigs, A);
    vec eigs_lcs = res.eigenvalues;
    // HydraPrint(eigs);
    // HydraPrint(eigs_lcs);
    REQUIRE(close(eigs, eigs_lcs)); // exact eigenvalue solution

    cx_mat tmatc(res.tmat, mat(N, N, fill::zeros));
    cx_mat A2 = cx_mat(V * tmatc * V.t());
    REQUIRE(norm(A - A2) < N * ortho_level);
  }

  {
    Log("  shastry sutherland N=16 ortho levels");

    std::string lfilename =
        HYDRA_DIRECTORY "/misc/data/shastry.16.HB.J.Jd.fsl.toml";
    auto lfile = FileToml(lfilename, 'r');
    auto bonds = BondList(lfile["Interactions"]);
    int n_sites = bonds.n_sites();
    bonds["J"] = 0.63;
    bonds["Jd"] = 1.00;
    auto group = PermutationGroup(lfile["Symmetries"]);

    std::vector<std::string> irrep_names = {"Gamma.C1.A", "M.C1.A", "X0.C1.A",
                                            "X1.C1.A"};

    for (int nup = 0; nup <= 6; ++nup) {
      for (std::string k : irrep_names) {
        auto irrep = Representation(lfile[k]);
        auto block = Spinhalf(n_sites, nup, group, irrep);

        Log("   nup: {} k: {}", nup, k);
        // HydraPrint(block);

        auto rstate = random_state(block, false);
        auto mult = [&bonds, &block](cx_vec const &v, cx_vec &w) {
          apply(bonds, block, v, block, w);
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

        // HydraPrint(res.num_reorthogonalizations);
        int n_iter2 = res.num_iterations;

        cx_mat V = res.V;
        // HydraPrint(res.num_iterations);
        // HydraPrint(cx_mat(V.t() * V ));
        double total_orthogonality = norm(V.t() * V - eye(n_iter2, n_iter2));
        if (n_iter > block.size()) {
          CHECK(total_orthogonality < 1e-6);
        } else {
          CHECK(total_orthogonality < n_iter * ortho_level);
        }
      }
    }
  }

  {
    Log("  shastry sutherland N=20 ortho levels");

    std::string lfilename =
        HYDRA_DIRECTORY "/misc/data/shastry.20.HB.J.Jd.fsl.toml";
    auto lfile = FileToml(lfilename, 'r');
    auto bonds = BondList(lfile["Interactions"]);
    int n_sites = bonds.n_sites();
    bonds["J"] = 0.63;
    bonds["Jd"] = 1.00;
    auto group = PermutationGroup(lfile["Symmetries"]);

    std::vector<std::string> irrep_names = {"Gamma.C1.A"};

    for (int nup = 0; nup <= 6; ++nup) {
      for (std::string k : irrep_names) {
        auto irrep = Representation(lfile[k]);
        auto block = Spinhalf(n_sites, nup, group, irrep);

        Log("   nup: {} k: {}", nup, k);
        // HydraPrint(block);

        auto rstate = random_state(block, false);
        auto mult = [&bonds, &block](cx_vec const &v, cx_vec &w) {
          apply(bonds, block, v, block, w);
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

        // HydraPrint(res.num_reorthogonalizations);
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
