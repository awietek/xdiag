#include "../../catch.hpp"

#include <iostream>
#include <xdiag/blocks/electron.hpp>
#include <xdiag/io/read.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;
using namespace xdiag::combinatorics;

// template <class bit_t>
// void test_indices(Electron<bit_t> const &electron) {
//   // int nsites = electron.nsites();

//   int64_t idx = 0;
//   for (bit_t ups : electron.reps_up_) {
//     auto const &dnss = electron.dns_for_up_rep(ups);
//     for (bit_t dns : dnss) {
//       // Log.out("{} {}", bits_to_string(ups, nsites),
//       // bits_to_string(dns, nsites));
//       int64_t idx2 = electron.index(ups, dns);
//       REQUIRE(idx == idx2);
//       ++idx;
//     }
//   }
// }

void test_electron_chain(int nsites) {
  using bits::bits_to_string;

  Log("electron_symmetric: chain block test. N: {}", nsites);

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
  auto group = PermutationGroup(permutation_array);

  // Create irrep with K momentum
  int64_t sum_of_dims = 0;

  for (int64_t nup = 0; nup <= nsites; ++nup) {
    for (int64_t ndn = 0; ndn <= nsites; ++ndn) {
      int64_t sum_of_dims_updn = 0;

      for (int64_t k = 0; k < nsites; ++k) {

        // Log.out("N: {}, nup: {}, ndn: {}, k:{} ", nsites, nup, ndn,
        // k);

        // Create characters
        std::vector<complex> chis;
        for (int64_t l = 0; l < nsites; ++l)
          chis.push_back({std::cos(2 * M_PI * l * k / nsites),
                          std::sin(2 * M_PI * l * k / nsites)});
        auto irrep = Representation(group, chis);
        auto electron2 = Electron(nsites, nup, ndn, irrep);

        sum_of_dims += electron2.size();
        sum_of_dims_updn += electron2.size();
      }
      CHECK(sum_of_dims_updn ==
            binomial(nsites, nup) * binomial(nsites, ndn));
    }
  }
  REQUIRE(sum_of_dims == pow(4, nsites));
}

TEST_CASE("electron_symmetric", "[electron]") {

  // Test the Hubbard chain
  for (int64_t nsites = 1; nsites < 7; ++nsites) {
    test_electron_chain(nsites);
  }

  // test a 3x3 triangular lattice
  int64_t nsites = 9;
  Log.out("electron_symmetric: triangular 3x3 block test");

  std::vector<std::pair<std::string, int64_t>> rep_name_mult = {
      {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
      {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
      {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
      {"Y.C1.A", 6}};

  std::string lfile =
      XDIAG_DIRECTORY "/misc/data/triangular.9.hop.sublattices.tsl.toml";
  auto fl = FileToml(lfile);

  int64_t sum_dim = 0;
  for (int64_t nup = 0; nup <= nsites; ++nup) {
    for (int64_t ndn = 0; ndn <= nsites; ++ndn) {
      int64_t sum_dim_updn = 0;

      for (auto [name, mult] : rep_name_mult) {
        auto irrep = read_representation(fl, name);
        auto electron = Electron(nsites, nup, ndn, irrep);
        int64_t dim = electron.size() * mult;
        // Log.out(
        //     "Hubbard Triangular 3x3: nsites: {}, nup: {}, ndn: {}, k: "
        //     "{}(x{}), size: {}",
        //     nsites, nup, ndn, name, mult, electron.size());
        // test_indices(electron);

        sum_dim_updn += dim;
        sum_dim += dim;
      }
      auto electron_nosym = Electron(nsites, nup, ndn);
      REQUIRE(sum_dim_updn == electron_nosym.size());
      REQUIRE(sum_dim_updn == binomial(nsites, nup) * binomial(nsites, ndn));

      // Log.out("size: {} {}", sum_dim_updn, electron_nosym.size());
    }
  }
  REQUIRE(sum_dim == (int64_t)pow(4, nsites));
}
