#include "../catch.hpp"

#include <iostream>

#include <hydra/all.h>

using namespace hydra;
using namespace hydra::combinatorics;

template <class bit_t> void print_states(Electron<bit_t> electron) {
  int n_sites = electron.n_sites();
  for (auto [ups, dns] : electron) {
    std::cout << bits_to_string(ups, n_sites) << ";"
              << bits_to_string(dns, n_sites) << "\n";
  }
}

template <class bit_t> void test_electron(int n_sites) {

  for (int nup = 0; nup <= n_sites; ++nup)
    for (int ndn = 0; ndn <= n_sites; ++ndn) {
      auto electron = Electron<bit_t>(n_sites, nup, ndn);
      // print_states(electron);
      // std::cout << "\n";
      bit_t cups = 0;
      bit_t cdns = 0;
      idx_t idx = 0;
      for (auto [ups, dns] : electron) {
        if (idx != 0) {
          REQUIRE(((dns > cdns) || (ups > cups)));
          REQUIRE(utils::popcnt(ups) == nup);
          REQUIRE(utils::popcnt(dns) == ndn);
        }
        cups = ups;
        cdns = dns;
        ++idx;
      }
      REQUIRE(electron.size() ==
              binomial(n_sites, nup) * binomial(n_sites, ndn));
      REQUIRE(electron.size() == idx);
    }

  // cyclic group as space group
  std::vector<std::vector<int>> permutations;
  for (int sym = 0; sym < n_sites; ++sym) {
    std::vector<int> permutation;
    for (int site = 0; site < n_sites; ++site) {
      int newsite = (site + sym) % n_sites;
      permutation.push_back(newsite);
    }
    permutations.push_back(permutation);
  }
  auto space_group = SpaceGroup<bit_t>(permutations);

  // HydraLog.out("n_sites: {}", n_sites);
  idx_t total_dim = 0;
  for (int nup = 0; nup <= n_sites; ++nup)
    for (int ndn = 0; ndn <= n_sites; ++ndn) {
      idx_t total_dim_nup_ndn = 0;

      for (int k = 0; k < n_sites; ++k) { // loop over momenta

        // Create irrep with K momentum
        std::vector<complex> chis;
        for (int l = 0; l < n_sites; ++l)
          chis.push_back({std::cos(2 * M_PI * l * k / n_sites),
                          std::sin(2 * M_PI * l * k / n_sites)});

        // HydraLog.out("k: {}", k);
        // for (auto chi : chis)
        //   std::cout << chi << " ";
        // std::cout << "\n";

        auto irrep = Representation(chis);
        auto electron = Electron<bit_t>(n_sites, nup, ndn, space_group, irrep);
        // print_states(electron);
        // std::cout << "\n";

        bit_t cups = 0;
        bit_t cdns = 0;
        idx_t idx = 0;
        for (auto [ups, dns] : electron) {
          if (idx != 0) {
            REQUIRE(((dns > cdns) || (ups > cups)));
            REQUIRE(utils::popcnt(ups) == nup);
            REQUIRE(utils::popcnt(dns) == ndn);
          }
          cups = ups;
          cdns = dns;
          ++idx;
        }
        REQUIRE(electron.size() == idx);
        // HydraLog.out("k: {}, nup: {}, ndn: {}, size: {}", k, nup, ndn, idx);
        total_dim_nup_ndn += idx;
        total_dim += idx;
      }
      idx_t dim = binomial(n_sites, nup) * binomial(n_sites, ndn);
      REQUIRE(total_dim_nup_ndn == dim);
      // HydraLog.out("totaldim_nup_ndn: {}, binomial: {}", total_dim_nup_ndn,
      //              dim);
    }
  idx_t dim = (idx_t)pow(4, n_sites);
  REQUIRE(total_dim == dim);
  // HydraLog.out("totaldim: {}, pow: {}\n", total_dim, dim);
}

TEST_CASE("electron", "[models]") {
  // test_electron<uint32>(4);

  for (int n_sites = 1; n_sites < 7; ++n_sites) {
    test_electron<hydra::uint16>(n_sites);
    test_electron<hydra::uint32>(n_sites);
    test_electron<hydra::uint64>(n_sites);
  }
}
