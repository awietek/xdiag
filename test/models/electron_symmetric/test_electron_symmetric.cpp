#include "../../catch.hpp"

#include <iostream>

#include <hydra/all.h>

using namespace hydra;
using namespace hydra::combinatorics;

template <class bit_t>
void test_indices(ElectronSymmetric<bit_t> const &electron) {
  for (auto [up, lower_upper] : electron.ups_lower_upper_) {
    idx_t lower = lower_upper.first;
    idx_t upper = lower_upper.second;
    for (idx_t idx = lower; idx < upper; ++idx) {
      bit_t dn = electron.dn(idx);

      auto [rep_up, rep_dn] = electron.representative(up, dn);
      REQUIRE(rep_up == up);
      REQUIRE(rep_dn == dn);
      idx_t idx2 = electron.index(up, dn);

      auto [rep_dn_switch, rep_up_switch] = electron.representative(dn, up);
      idx_t idx_switch = electron.index_switch(rep_up_switch, rep_dn_switch);
      idx_t idx3 = electron.index_switch_to_index(idx_switch);

      REQUIRE(idx == idx2);
      REQUIRE(idx == idx3);
    }
  }
}

template <class bit_t>
void test_representative_character(ElectronSymmetric<bit_t> const &electron) {
  int n_sites = electron.n_sites();
  int nup = electron.n_up();
  int ndn = electron.n_dn();
  auto space_group = electron.group_action();
  auto irrep = electron.irrep();
  for (auto up : Combinations(n_sites, nup)) {
    for (auto dn : Combinations(n_sites, ndn)) {

      auto [rep_up, rep_dn, sym] = electron.representative_index(up, dn);
      idx_t idx1 = electron.index(rep_up, rep_dn);

      auto [rep_dn_switch, rep_up_switch, sym_switch] =
          electron.representative_index(dn, up);
      idx_t idx_switch = electron.index_switch(rep_up_switch, rep_dn_switch);
      idx_t idx2 = electron.index_switch_to_index(idx_switch);
      REQUIRE(idx1 == idx2);

      if (index_valid(idx1)) {
        complex chi1 = irrep.character(sym) * space_group.fermi_sign(sym, up) *
                       space_group.fermi_sign(sym, dn);
        complex chi2 = irrep.character(sym_switch) *
                       space_group.fermi_sign(sym_switch, up) *
                       space_group.fermi_sign(sym_switch, dn) *
                       electron.character_switch_[idx_switch];
        REQUIRE(lila::close(chi1, chi2));
      }
    }
  }
}

template <class bit_t> void test_electron_chain(int n_sites) {

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
  auto space_group = PermutationGroup(permutations);

  // Create irrep with K momentum
  idx_t sum_of_dims = 0;

  for (int nup = 0; nup <= n_sites; ++nup) {
    for (int ndn = 0; ndn <= n_sites; ++ndn) {
      idx_t sum_of_dims_updn = 0;

      for (int k = 0; k < n_sites; ++k) {

        lila::Log.out("N: {}, nup: {}, ndn: {}, k:{} ", n_sites, nup, ndn, k);

        // Create characters
        std::vector<complex> chis;
        for (int l = 0; l < n_sites; ++l)
          chis.push_back({std::cos(2 * M_PI * l * k / n_sites),
                          std::sin(2 * M_PI * l * k / n_sites)});
        auto irrep = Representation(chis);

        auto electron1 =
            ElectronSymmetric<bit_t>(n_sites, nup, ndn, space_group, irrep);
        auto electron2 =
            ElectronSymmetricV2<bit_t>(n_sites, nup, ndn, space_group, irrep);

        sum_of_dims += electron2.size();
        sum_of_dims_updn += electron2.size();
	

	// // Compare up, dn order
        // {
        //   std::vector<bit_t> ups1;
        //   std::vector<bit_t> dns1;
        //   std::vector<double> norms1;
        //   for (auto [up, lower_upper] : electron1.ups_lower_upper_) {
        //     idx_t lower = lower_upper.first;
        //     idx_t upper = lower_upper.second;

        //     for (idx_t idx = lower; idx < upper; ++idx) {
        //       bit_t dn = electron1.dns_[idx];
        //       double norm = electron1.norms_[idx];
        //       ups1.push_back(up);
        //       dns1.push_back(dn);
        //       norms1.push_back(norm);
        //     }
        //   }

        //   std::vector<bit_t> ups2;
        //   std::vector<bit_t> dns2;
        //   std::vector<double> norms2;
        //   for (bit_t ups : electron2.reps_up_) {

        //     auto const &dnss = electron2.dns_for_up_rep(ups);
        //     auto const &norms = electron2.norms_for_up_rep(ups);

        //     for (idx_t idx = 0; idx < dnss.size(); ++idx) {
        //       bit_t dns = dnss[idx];
        //       double norm = norms[idx];

        //       ups2.push_back(ups);
        //       dns2.push_back(dns);
        //       norms2.push_back(norm);
        //     }
        //   }

        //   REQUIRE(ups1.size() == dns1.size());
        //   REQUIRE(ups1.size() == norms1.size());
        //   REQUIRE(ups1.size() == ups2.size());
        //   REQUIRE(ups1.size() == dns2.size());
        //   REQUIRE(ups1.size() == norms2.size());

        //   for (idx_t idx = 0; idx < ups1.size(); ++idx) {
        //     lila::Log.out("{};{} {:.4f}  <->  {};{} {:.4f}",
        //                   bits_to_string(ups1[idx], n_sites),
        //                   bits_to_string(dns1[idx], n_sites), norms1[idx],
        //                   bits_to_string(ups2[idx], n_sites),
        //                   bits_to_string(dns2[idx], n_sites), norms2[idx]);
        //     REQUIRE(ups1[idx] == ups2[idx]);
        //     REQUIRE(dns1[idx] == dns2[idx]);
        //     REQUIRE(norms1[idx] == norms2[idx]);
        //   }
        // }

	// Compare dn, up order
        {
          std::vector<bit_t> ups1;
          std::vector<bit_t> dns1;
          std::vector<double> norms1;
          for (auto [dn, lower_upper] : electron1.dns_lower_upper_) {
            idx_t lower = lower_upper.first;
            idx_t upper = lower_upper.second;

            for (idx_t idx = lower; idx < upper; ++idx) {
              bit_t up = electron1.ups_[idx];
              double norm = electron1.norms_[electron1.index_switch_to_index(idx)];
              ups1.push_back(up);
              dns1.push_back(dn);
              norms1.push_back(norm);
            }
          }

          std::vector<bit_t> ups2;
          std::vector<bit_t> dns2;
          std::vector<double> norms2;
          for (bit_t dns : electron2.reps_dn_) {

            auto const &upss = electron2.ups_for_dn_rep(dns);
            auto const &norms = electron2.norms_for_dn_rep(dns);

            for (idx_t idx = 0; idx < upss.size(); ++idx) {
              bit_t ups = upss[idx];
              double norm = norms[idx];

              ups2.push_back(ups);
              dns2.push_back(dns);
              norms2.push_back(norm);
            }
          }

          REQUIRE(ups1.size() == dns1.size());
          REQUIRE(ups1.size() == norms1.size());
          REQUIRE(ups1.size() == ups2.size());
          REQUIRE(ups1.size() == dns2.size());
          REQUIRE(ups1.size() == norms2.size());

          for (idx_t idx = 0; idx < ups1.size(); ++idx) {
            lila::Log.out("{};{} {:.4f}  <->  {};{} {:.4f}",
                          bits_to_string(ups1[idx], n_sites),
                          bits_to_string(dns1[idx], n_sites), norms1[idx],
                          bits_to_string(ups2[idx], n_sites),
                          bits_to_string(dns2[idx], n_sites), norms2[idx]);
            REQUIRE(ups1[idx] == ups2[idx]);
            REQUIRE(dns1[idx] == dns2[idx]);
            REQUIRE(norms1[idx] == norms2[idx]);
          }
        }






        // if (electron.size() > 0) {
        //   // test_indices(electron);
        //   // test_representative_character(electron);
        // }
      }
      CHECK(sum_of_dims_updn ==
            binomial(n_sites, nup) * binomial(n_sites, ndn));
    }
  }
  REQUIRE(sum_of_dims == pow(4, n_sites));
}

TEST_CASE("ElectronSymmetric", "[models]") {

  // Test the Hubbard chain
  for (int n_sites = 1; n_sites < 7; ++n_sites) {
    test_electron_chain<hydra::uint16>(n_sites);
    // test_electron_chain<hydra::uint32>(n_sites);
    // test_electron_chain<hydra::uint64>(n_sites);
  }

  // // test a 3x3 triangular lattice
  // int n_sites = 9;

  // std::vector<std::pair<std::string, int>> rep_name_mult = {
  //     {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
  //     {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
  //     {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
  //     {"Y.C1.A", 6}};

  // std::string lfile = "data/triangular.9.hop.sublattices.tsl.lat";
  // auto permutations = hydra::read_permutations(lfile);
  // auto space_group = PermutationGroup(permutations);

  // idx_t sum_dim = 0;
  // for (int nup = 0; nup <= n_sites; ++nup) {
  //   for (int ndn = 0; ndn <= n_sites; ++ndn) {
  //     idx_t sum_dim_updn = 0;

  //     for (auto [name, mult] : rep_name_mult) {
  //       auto irrep = read_represenation(lfile, name);
  //       auto electron =
  //           ElectronSymmetric<uint16>(n_sites, nup, ndn, space_group,
  //           irrep);

  //       idx_t dim = electron.size() * mult;
  //       // lila::Log.out(
  //       //     "Hubbard Triangular 3x3: n_sites: {}, nup: {}, ndn: {}, k: "
  //       //     "{}(x{}), size: {}",
  //       //     n_sites, nup, ndn, name, mult, electron.size());
  //       test_indices(electron);
  //       test_representative_character(electron);

  //       sum_dim_updn += dim;
  //       sum_dim += dim;
  //     }
  //     auto electron_nosym = Electron<uint16>(n_sites, nup, ndn);
  //     REQUIRE(sum_dim_updn == electron_nosym.size());
  //     REQUIRE(sum_dim_updn == binomial(n_sites, nup) * binomial(n_sites,
  //     ndn));

  //     // lila::Log.out("size: {} {}", sum_dim_updn, electron_nosym.size());
  //   }
  // }
  // REQUIRE(sum_dim == (idx_t)pow(4, n_sites));

  // int n_sites = 6;
  // int nup = 3;
  // int ndn = 3;

  // // cyclic group as space group
  // std::vector<std::vector<int>> permutations;
  // for (int sym = 0; sym < n_sites; ++sym) {
  //   std::vector<int> permutation;
  //   for (int site = 0; site < n_sites; ++site) {
  //     int newsite = (site + sym) % n_sites;
  //     permutation.push_back(newsite);
  //   }
  //   permutations.push_back(permutation);
  // }
  // auto space_group = PermutationGroup(permutations);

  // int k = 0;

  // // Create characters
  // std::vector<complex> chis;
  // for (int l = 0; l < n_sites; ++l)
  //   chis.push_back({std::cos(2 * M_PI * l * k / n_sites),
  //                   std::sin(2 * M_PI * l * k / n_sites)});
  // auto irrep = Representation(chis);

  // ElectronSymmetricV2(n_sites, nup, ndn, space_group, irrep);
}
