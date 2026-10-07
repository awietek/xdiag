// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cstring>

#include <tests/blocks/electron/testcases_electron.hpp>
#include <tests/blocks/tj/testcases_tj.hpp>
#include <tests/catch.hpp>

#include <xdiag/linalg/lanczos/eigs_lanczos.hpp>
#include <xdiag/states/apply.hpp>
#include <xdiag/states/fill.hpp>
#include <xdiag/states/random_state.hpp>
#include <xdiag/utils/logger.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace xdiag;

// eigs_lanczos runs the Lanczos recurrence twice: the first pass yields the
// tridiagonal matrix, the second pass regenerates the very same Krylov vectors
// and accumulates the eigenvector from them. This is only correct if every
// matrix-vector product and inner product is bitwise reproducible. Far beyond
// convergence the recurrence amplifies any roundoff difference exponentially
// (loss of orthogonality, ghost eigenvalues), so these tests deliberately run
// many more iterations than needed for convergence. Cf. GitHub discussion #119.

template <typename coeff_t>
static bool bitwise_equal(arma::Col<coeff_t> const &a,
                          arma::Col<coeff_t> const &b) {
  if (a.n_elem != b.n_elem) {
    return false;
  }
  return std::memcmp(a.memptr(), b.memptr(), a.n_elem * sizeof(coeff_t)) == 0;
}

static void require_bitwise_equal(EigsLanczosResult const &r1,
                                  EigsLanczosResult const &r2) {
  REQUIRE(r1.niterations == r2.niterations);
  REQUIRE(bitwise_equal(r1.alphas, r2.alphas));
  REQUIRE(bitwise_equal(r1.betas, r2.betas));
  REQUIRE(bitwise_equal(r1.eigenvalues, r2.eigenvalues));
  if (isreal(r1.eigenvectors)) {
    REQUIRE(bitwise_equal(r1.eigenvectors.vector(0, false),
                          r2.eigenvectors.vector(0, false)));
  } else {
    REQUIRE(bitwise_equal(r1.eigenvectors.vectorC(0, false),
                          r2.eigenvectors.vectorC(0, false)));
  }
}

// The reconstructed vector must be an eigenvector of ops to the eigenvalue
// found in the first pass: check the energy and the residual after
// normalization.
static void require_eigenvector(OpSum const &ops, EigsLanczosResult const &r) {
  double e0 = r.eigenvalues(0);
  State v = r.eigenvectors.col(0);
  State Hv = v;
  apply(ops, v, Hv);
  if (isreal(v)) {
    arma::vec vv = v.vector(0, false);
    arma::vec hv = Hv.vector(0, false);
    double nrm = arma::norm(vv);
    REQUIRE(nrm > 1e-8);
    double e = arma::dot(vv, hv) / (nrm * nrm);
    double residual = arma::norm(hv - e0 * vv) / nrm;
    REQUIRE(std::abs(e - e0) < 1e-8);
    REQUIRE(residual < 1e-8);
  } else {
    arma::cx_vec vv = v.vectorC(0, false);
    arma::cx_vec hv = Hv.vectorC(0, false);
    double nrm = arma::norm(vv);
    REQUIRE(nrm > 1e-8);
    double e = std::real(arma::cdot(vv, hv)) / (nrm * nrm);
    double residual = arma::norm(hv - e0 * vv) / nrm;
    REQUIRE(std::abs(e - e0) < 1e-8);
    REQUIRE(residual < 1e-8);
  }
}

static void check_rerun(OpSum const &ops, Block const &block,
                        int64_t niterations) {
  // precision 0: the convergence criterion is never met, so exactly
  // niterations Lanczos steps are performed in both passes.
  double precision = 0.0;

  // Fixed start vector, so the runs below differ in nothing but roundoff.
  State state0(block, isreal(ops) && isreal(block));
  fill(state0, RandomState(42));

  Log("  rerun with identical settings");
  EigsLanczosResult r1 = eigs_lanczos(ops, state0, 1, precision, niterations);
  EigsLanczosResult r2 = eigs_lanczos(ops, state0, 1, precision, niterations);
  REQUIRE(r1.niterations == niterations);
  require_bitwise_equal(r1, r2);
  require_eigenvector(ops, r1);

#ifdef _OPENMP
  // The result must not depend on the number of threads either.
  int nthreads_before = omp_get_max_threads();
  int nthreads = std::min(omp_get_num_procs(), 16);
  if (nthreads > 1) {
    Log("  rerun with 1 and {} threads", nthreads);
    omp_set_num_threads(1);
    EigsLanczosResult rs = eigs_lanczos(ops, state0, 1, precision, niterations);
    omp_set_num_threads(nthreads);
    EigsLanczosResult rp = eigs_lanczos(ops, state0, 1, precision, niterations);
    require_bitwise_equal(rs, rp);
    require_bitwise_equal(r1, rp);
  }
  omp_set_num_threads(nthreads_before);
#endif
}

TEST_CASE("eigs_lanczos_rerun_bitwise", "[lanczos]") {
  using namespace xdiag::testcases::tj;
  using xdiag::testcases::electron::get_cyclic_group_irreps;

  // Without lattice symmetries every output element receives at most one
  // contribution per term, so the order of the atomic updates in the scatter
  // kernels never mattered. Still a useful guarantee.
  {
    int nsites = 10;
    int nup = 4;
    int ndn = 3;
    int64_t niterations = 300; // far beyond convergence (dim = 4200)

    Log("eigs_lanczos_rerun_bitwise: real tJ chain");
    OpSum ops = tJchain(nsites, 1.0, 0.3);
    tJ block(nsites, nup, ndn);
    check_rerun(ops, block, niterations);
  }

  // With lattice symmetries several representatives map onto the same output
  // representative within one term, so several threads add to the same
  // element and the scatter kernels are order dependent. This is the setting
  // of GitHub discussion #119 (16-site t-J chain, 8 up, 6 dn, momentum sector,
  // dimension 22,470, 500 Lanczos steps).
  {
    int nsites = 16;
    int nup = 8;
    int ndn = 6;
    int64_t niterations = 500;
    OpSum ops = tJchain(nsites, 1.0, 0.3);
    std::vector<Representation> irreps = get_cyclic_group_irreps(nsites);

    Log("eigs_lanczos_rerun_bitwise: symmetric tJ chain, k = 0 (real)");
    {
      tJ block(nsites, nup, ndn, irreps[0]);
      check_rerun(ops, block, niterations);
    }

    Log("eigs_lanczos_rerun_bitwise: symmetric tJ chain, k = 1 (complex)");
    {
      tJ block(nsites, nup, ndn, irreps[1]);
      check_rerun(ops, block, niterations);
    }
  }
}
