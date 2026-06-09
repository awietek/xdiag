// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../catch.hpp"

#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <xdiag/algebra/symmetrize.hpp>
#include <xdiag/armadillo.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/boson.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/config.hpp>
#include <xdiag/io/file_toml.hpp>
#include <xdiag/io/read.hpp>
#include <xdiag/linalg/sparse_diag.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/math/isapprox.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/correlation_matrix.hpp>
#include <xdiag/states/inner.hpp>
#include <xdiag/states/state.hpp>
#include <xdiag/symmetries/cyclic_group.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/representation.hpp>
#include <xdiag/utils/error.hpp>

using namespace xdiag;

// Naive correlation matrix: C(i,j) = <psi| build(i,j) |psi> for an OpSum
// builder. On a symmetry-adapted block, build(i,j) must already be symmetrized
// so the bare inner / apply stays inside the block.
template <typename Builder>
static arma::mat naive_correlation(State const &psi, Builder build) {
  int64_t n = psi.nsites();
  arma::mat r(n, n);
  for (int64_t i = 0; i < n; ++i) {
    for (int64_t j = 0; j < n; ++j) {
      r(i, j) = inner(psi, build(i, j), psi);
    }
  }
  return r;
}

// Complex variant (innerC) for complex states.
template <typename Builder>
static arma::cx_mat naive_correlationC(State const &psi, Builder build) {
  int64_t n = psi.nsites();
  arma::cx_mat r(n, n);
  for (int64_t i = 0; i < n; ++i) {
    for (int64_t j = 0; j < n; ++j) {
      r(i, j) = innerC(psi, build(i, j), psi);
    }
  }
  return r;
}

// Kitaev-Gamma model on the N=8 honeycomb: a genuinely complex Hamiltonian
// (complex Matrix ops on bond-dependent X/Y/Z terms). Loaded from the shared
// lattice file, mirroring tests/blocks/spinhalf/test_spinhalf_strategies.
static OpSum kitaev_gamma_opsum(double K, double G) {
  std::string lfile = XDIAG_DIRECTORY
      "/misc/data/kitaev_gamma/lattice-files/"
      "honeycomb.8.HeisenbergKitaevGamma.fsl.toml";
  auto fl = FileToml(lfile);
  auto ops_read = read_opsum(fl, "Interactions");

  arma::mat zero2(2, 2, arma::fill::zeros);
  arma::cx_mat sx(arma::mat({{0., 0.5}, {0.5, 0.}}), zero2);
  arma::cx_mat sy(zero2, arma::mat({{0., -0.5}, {0.5, 0.}}));
  arma::cx_mat sz(arma::mat({{0.5, 0.0}, {0.0, -0.5}}), zero2);
  arma::cx_mat sxsx = arma::kron(sx, sx);
  arma::cx_mat sysy = arma::kron(sy, sy);
  arma::cx_mat szsz = arma::kron(sz, sz);
  arma::cx_mat gsx = arma::kron(sy, sz) + arma::kron(sz, sy);
  arma::cx_mat gsy = arma::kron(sx, sz) + arma::kron(sz, sx);
  arma::cx_mat gsz = arma::kron(sx, sy) + arma::kron(sy, sx);

  OpSum ops;
  for (auto term : ops_read) {
    Op const &op = term.monomial[0];
    std::string type = op.type();
    auto sites = op.sites();
    if (type == "KITAEVX") {
      ops += K * Op("Matrix", sites, sxsx);
    } else if (type == "KITAEVY") {
      ops += K * Op("Matrix", sites, sysy);
    } else if (type == "KITAEVZ") {
      ops += K * Op("Matrix", sites, szsz);
    } else if (type == "GAMMAX") {
      ops += G * Op("Matrix", sites, gsx);
    } else if (type == "GAMMAY") {
      ops += G * Op("Matrix", sites, gsy);
    } else if (type == "GAMMAZ") {
      ops += G * Op("Matrix", sites, gsz);
    }
  }
  return ops;
}

TEST_CASE("correlation_matrix", "[states]") try {
  // ===== Spin-1/2 Heisenberg ring: <S+_i S-_j> =====
  {
    int64_t n = 4;
    OpSum H;
    for (int64_t i = 0; i < n; ++i) {
      H += Op("SdotS", {i, (i + 1) % n});
    }

    // No permutation symmetry (nup sector only): bare ops act in-block.
    {
      auto block = Spinhalf(n, n / 2);
      State psi = std::get<1>(eig0(H, block));
      arma::mat c = correlation_matrix(psi, "S+", "S-");
      arma::mat naive = naive_correlation(psi, [&](int64_t i, int64_t j) {
        return OpSum(Monomial{Op("S+", i), Op("S-", j)});
      });
      REQUIRE(isapprox(c, naive, 1e-10, 1e-10));
    }

    // Translation symmetry (k = 0): correlation_matrix symmetrizes internally;
    // the naive uses the same group-symmetrized operator.
    {
      auto group = cyclic_group(n);
      auto block = Spinhalf(n, n / 2, cyclic_group_irrep(n, 0));
      State psi = std::get<1>(eig0(H, block));
      arma::mat c = correlation_matrix(psi, "S+", "S-");
      arma::mat naive = naive_correlation(psi, [&](int64_t i, int64_t j) {
        return symmetrize(OpSum(Monomial{Op("S+", i), Op("S-", j)}), group);
      });
      REQUIRE(isapprox(c, naive, 1e-10, 1e-10));
    }
  }

  // ===== Bose-Hubbard ring: <Adag_i A_j> (single-particle density matrix) =====
  {
    int64_t n = 4, d = 3, number = 3;
    OpSum H;
    for (int64_t i = 0; i < n; ++i) {
      H += Op("Hop", {i, (i + 1) % n});
    }
    H += 2.0 * Op("HubbardU");

    {
      auto block = Boson(n, d, number);
      State psi = std::get<1>(eig0(H, block));
      arma::mat c = correlation_matrix(psi, "Adag", "A");
      arma::mat naive = naive_correlation(psi, [&](int64_t i, int64_t j) {
        return OpSum(Monomial{Op("Adag", i), Op("A", j)});
      });
      REQUIRE(isapprox(c, naive, 1e-10, 1e-10));
    }

    {
      auto group = cyclic_group(n);
      auto block = Boson(n, d, number, cyclic_group_irrep(n, 0));
      State psi = std::get<1>(eig0(H, block));
      arma::mat c = correlation_matrix(psi, "Adag", "A");
      arma::mat naive = naive_correlation(psi, [&](int64_t i, int64_t j) {
        return symmetrize(OpSum(Monomial{Op("Adag", i), Op("A", j)}), group);
      });
      REQUIRE(isapprox(c, naive, 1e-10, 1e-10));
    }
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
}

// Complex model: Kitaev-Gamma honeycomb. Exercises the complex paths
// (correlation_matrixC / dotC), the real-ness guard, and -- via the
// non-Hermitian S+/S- correlator -- the hermitian-conjugate dedup branch
// (op(j,i) = op(i,j)^dagger, value conjugated). Both with and without a C2
// point-group symmetry, validated against the naive innerC loop.
TEST_CASE("correlation_matrix_kitaev", "[states]") try {
  double K = -.30901699437494742412;
  double G = .95105651629515357210;
  OpSum ops = kitaev_gamma_opsum(K, G);

  // Hermitian correlators (Sx,Sx),(Sy,Sy),(Sz,Sz) and the non-Hermitian
  // (S+,S-) which drives the conjugate dedup branch.
  std::vector<std::pair<std::string, std::string>> pairs = {
      {"Sx", "Sx"}, {"Sy", "Sy"}, {"Sz", "Sz"}, {"S+", "S-"}};

  // --- without permutation symmetry ---
  {
    auto block = Spinhalf(8);
    State psi = std::get<1>(eig0(ops, block));

    // a real correlation matrix cannot be formed from a complex state
    REQUIRE_THROWS(correlation_matrix(psi, "Sz", "Sz"));

    for (auto const &[t1, t2] : pairs) {
      arma::cx_mat c = correlation_matrixC(psi, t1, t2);
      arma::cx_mat naive = naive_correlationC(psi, [&](int64_t i, int64_t j) {
        return OpSum(Monomial{Op(t1, i), Op(t2, j)});
      });
      REQUIRE(isapprox(c, naive, 1e-8, 1e-8));
    }
  }

  // --- with a C2 point-group symmetry ---
  {
    std::string lfile = XDIAG_DIRECTORY
        "/misc/data/kitaev_gamma/lattice-files/"
        "honeycomb.8.HeisenbergKitaevGamma.fsl.toml";
    auto fl = FileToml(lfile);
    auto irrep = read_representation(fl, "Gamma.C2.A");
    auto group = irrep.group();
    auto block = Spinhalf(8, irrep);
    State psi = std::get<1>(eig0(ops, block));

    for (auto const &[t1, t2] : pairs) {
      arma::cx_mat c = correlation_matrixC(psi, t1, t2);
      arma::cx_mat naive = naive_correlationC(psi, [&](int64_t i, int64_t j) {
        return symmetrize(OpSum(Monomial{Op(t1, i), Op(t2, j)}), group);
      });
      REQUIRE(isapprox(c, naive, 1e-8, 1e-8));
    }
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
}
