// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <iostream>

#include <tests/catch.hpp>

#include <xdiag/bits/popcount.hpp>
#include <xdiag/bits/zero_one.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/symmetries/fermi_sign.hpp>
#include <xdiag/symmetries/permutation.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;

// Slow but simple reference implementation we test against
template <class bit_t>
static bool fermi_bool_of_permutation_sort(bit_t state,
                                           Permutation const &permutation) {
  int nsites = permutation.size();
  std::vector<int64_t> work(nsites * 2);

  int64_t *iota = work.data();
  int64_t *to = work.data() + bits::popcount(state);

  // find out where fermions are mapped to
  int64_t n_fermion = 0;
  for (int64_t site = 0; bits::nonzero(state); ++site) {
    if (bits::nonzero(state & bits::one<bit_t>(nsites))) {
      iota[n_fermion] = n_fermion;
      to[n_fermion++] = permutation[site];
    }
    state >>= 1;
  }

  // Find sorting permutation -> iota
  std::sort(iota, iota + n_fermion, [&to](const int64_t &a, const int64_t &b) {
    return to[a] < to[b];
  });

  // compute sign in O(n_fermions) by cycle decomposition
  bool sign = false;
  bit_t visited =
      (bits::one<bit_t>(nsites) << n_fermion) - bits::one<bit_t>(nsites);
  int64_t next, L;
  bit_t mask;
  for (int64_t k = 0; k < n_fermion; ++k) {
    if (bits::nonzero(visited & (bits::one<bit_t>(nsites) << k))) {
      next = k;
      L = 0;
      mask = (bits::one<bit_t>(nsites) << next);
      while (bits::nonzero(visited & mask)) {
        ++L;
        visited ^= mask;
        next = iota[next];
        mask = (bits::one<bit_t>(nsites) << next);
      }
      if (!(L & 1))
        sign = !sign;
    }
  }
  return sign;
}

template <class bit_t> void test_fermi_sign_subsets(int64_t nsites) {
  using combinatorics::Subsets;

  // test cyclic group
  std::vector<Permutation> permutation_array;
  for (int64_t sym = 0; sym < nsites; ++sym) {

    std::vector<int64_t> pv;
    for (int64_t site = 0; site < nsites; ++site) {
      int64_t newsite = (site + sym) % nsites;
      pv.push_back(newsite);
    }
    permutation_array.push_back(Permutation(pv));
  }

  // auto sym_op = SpaceGroupOperator<bit_t>(nsites, permutation_array);

  for (auto state : Subsets<bit_t>(nsites)) {
    for (int64_t sym = 0; sym < nsites; ++sym) {

      // auto tstate = sym_op.apply(sym, state);
      // std::cout << "sym " << sym << "\n";
      // std::cout << bits_to_string(state, nsites) << "\n";
      // std::cout << bits_to_string(tstate, nsites) << "\n";

      auto fermi1 =
          symmetries::fermi_bool_of_permutation(state, permutation_array[sym]);
      auto fermi2 =
          fermi_bool_of_permutation_sort(state, permutation_array[sym]);

      // std::cout << fermi1 << " " << fermi2 << "\n";
      // std::cout << "\n";
      REQUIRE(fermi1 == fermi2);
    }
    // std::cout << "-----\n";
  }
}

template <class bit_t> void test_fermi_sign_combinations(int64_t nsites) {
  using combinatorics::Combinations;

  // test cyclic group
  std::vector<Permutation> permutation_array;
  for (int64_t sym = 0; sym < nsites; ++sym) {

    std::vector<int64_t> pv;
    for (int64_t site = 0; site < nsites; ++site) {
      int64_t newsite = (site + sym) % nsites;
      pv.push_back(newsite);
    }
    permutation_array.push_back(Permutation(pv));
  }

  for (int k = 0; k < nsites; ++k) {
    for (auto state : Combinations<bit_t>(nsites, k)) {
      for (int64_t sym = 0; sym < nsites; ++sym) {

        auto fermi1 = symmetries::fermi_bool_of_permutation(
            state, permutation_array[sym]);
        auto fermi2 =
            fermi_bool_of_permutation_sort(state, permutation_array[sym]);
        REQUIRE(fermi1 == fermi2);
      }
      // std::cout << "-----\n";
    }
  }
}

TEST_CASE("fermi_sign", "[symmetries]") {
  using namespace xdiag;
  using namespace bits;

  Log("Test fermi_sign");
  for (int64_t nsites = 1; nsites < 8; ++nsites) {
    test_fermi_sign_subsets<uint32_t>(nsites);
    test_fermi_sign_subsets<uint64_t>(nsites);
    test_fermi_sign_combinations<BitsetDynamic>(nsites);
    test_fermi_sign_combinations<BitsetStatic1>(nsites);
    test_fermi_sign_combinations<BitsetStatic2>(nsites);
    test_fermi_sign_combinations<BitsetStatic4>(nsites);
    test_fermi_sign_combinations<BitsetStatic8>(nsites);
  }
}
