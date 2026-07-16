// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <string>
#include <vector>

#include <xdiag/basis/basis_symmetric.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/config.hpp>
#include <xdiag/io/file_toml.hpp>
#include <xdiag/io/read.hpp>
#include <xdiag/symmetries/action/site_permutation.hpp>
#include <xdiag/symmetries/cyclic_group.hpp>
#include <xdiag/symmetries/representation.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;

// Checks representative_data consistency for a single state
template <typename bit_t, class Basis>
void check_state(bit_t state, Basis const &basis) {
  auto action = symmetries::SitePermutation(basis.group());
  auto [raw, sym, norm_out] = basis.representative_data(state);

  if (raw == 0)
    return;

  int64_t idx = raw - 1;
  bit_t rep = basis[idx];

  REQUIRE(action.apply(sym, state) == rep);
  REQUIRE(norm_out == Approx(basis.norm(idx)));
}

// Checks that the iterator visits all representatives exactly once, in order,
// and that each representative maps to itself under representative_data
template <class Basis> void check_iterator(Basis const &basis) {
  int64_t count = 0;
  for (auto state : basis) {
    auto [raw, sym, norm_out] = basis.representative_data(state);
    REQUIRE(raw != 0);
    REQUIRE(basis[raw - 1] == state);
    REQUIRE(raw - 1 == count);
    ++count;
  }
  REQUIRE(count == basis.size());
}

template <typename bit_t> void test_no_sz(Representation const &irrep) {
  using namespace basis;
  using namespace combinatorics;

  int64_t nsites = irrep.group().nsites();
  auto enumeration = Subsets<bit_t>(nsites);
  auto b = BasisSymmetric<Subsets<bit_t>>(enumeration, irrep.group(),
                                          irrep.characters());

  for (auto state : enumeration)
    check_state(state, b);
  check_iterator(b);
}

template <typename bit_t>
void test_sz(int64_t nup, Representation const &irrep) {
  using namespace basis;
  using namespace combinatorics;

  int64_t nsites = irrep.group().nsites();
  auto enumeration = Combinations<bit_t>(nsites, nup);
  auto b = BasisSymmetric<Combinations<bit_t>>(enumeration, irrep.group(),
                                               irrep.characters());

  for (auto state : enumeration)
    check_state(state, b);
  check_iterator(b);
}

template <typename bit_t> void test_basis_symmetric_cyclic() {
  for (int64_t nsites = 1; nsites <= 6; ++nsites) {
    for (int64_t k = 0; k < nsites; ++k) {
      auto irrep = cyclic_group_irrep(nsites, k);
      test_no_sz<bit_t>(irrep);
      for (int64_t nup = 0; nup <= nsites; ++nup)
        test_sz<bit_t>(nup, irrep);
    }
  }
}

template <typename bit_t> void test_basis_symmetric_lattice() {
  {
    Log("basis_symmetric: triangular 3x3");
    int64_t nsites = 9;
    std::string lfile = XDIAG_DIRECTORY
        "/misc/data/triangular.9.Jz1Jz2Jx1Jx2D1.sublattices.tsl.toml";
    auto fl = FileToml(lfile);
    std::vector<std::string> irrep_names = {
        "Gamma.D6.A1", "Gamma.D6.A2", "Gamma.D6.B1", "Gamma.D6.B2",
        "Gamma.D6.E1", "Gamma.D6.E2", "K.D3.A1",     "K.D3.A2",
        "K.D3.E",      "Y.D1.A",      "Y.D1.B"};
    for (auto const &name : irrep_names) {
      auto irrep = read_representation(fl, name);
      test_no_sz<bit_t>(irrep);
      for (int64_t nup = 0; nup <= nsites; ++nup)
        test_sz<bit_t>(nup, irrep);
    }
  }

  {
    Log("basis_symmetric: square 2x4");
    int64_t nsites = 8;
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/square.8.heisenberg.2sl.toml";
    auto fl = FileToml(lfile);
    std::vector<std::string> irrep_names = {
        "Gamma.D4.A1", "Gamma.D4.A2", "Gamma.D4.B1", "Gamma.D4.B2",
        "Gamma.D4.E",  "M.D4.A1",     "M.D4.A2",     "M.D4.B1",
        "M.D4.B2",     "M.D4.E",      "Sigma.D1.A",  "Sigma.D1.B",
        "X.D2.A1",     "X.D2.A2",     "X.D2.B1",     "X.D2.B2"};
    for (auto const &name : irrep_names) {
      auto irrep = read_representation(fl, name);
      test_no_sz<bit_t>(irrep);
      for (int64_t nup = 0; nup <= nsites; ++nup)
        test_sz<bit_t>(nup, irrep);
    }
  }
}

template <typename bit_t> void test_basis_symmetric() {
  Log("cyclic group irreps");
  test_basis_symmetric_cyclic<bit_t>();
  Log("lattice data");
  test_basis_symmetric_lattice<bit_t>();
}

TEST_CASE("basis_symmetric", "[basis]") try {
  Log("Test basis_symmetric");
  Log("uint32_t");
  test_basis_symmetric<uint32_t>();
  Log("uint64_t");
  test_basis_symmetric<uint64_t>();
  Log("Done");
} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}
