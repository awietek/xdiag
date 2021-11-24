#include "../../catch.hpp"

#include <iostream>

#include "testcases_tj.h"
#include <hydra/all.h>
#include <lila/all.h>

using namespace hydra;

TEST_CASE("tj_utils", "[blocks][tj]") {
  using namespace hydra::combinatorics;
  using namespace hydra::utils;
  using namespace hydra::bitops;
  lila::Log("tj_utils: test");

  for (int N = 0; N <= 6; ++N)
    for (int nup = 0; nup <= N; ++nup)
      for (int ndn = 0; ndn <= N - nup; ++ndn) {
        // lila::Log.out("N: {},  nup: {}, ndn: {}", N, nup, ndn);

        for (auto up : Combinations(N, nup))
          for (auto dnc : Combinations(N - nup, ndn)) {
            auto [upc, dn] = up_dnc_to_upc_dn(up, dnc);
            auto [up2, dnc2] = upc_dn_to_up_dnc(upc, dn);
            REQUIRE(up == up2);
            REQUIRE(dnc == dnc2);
            // lila::Log.out("up : {} dnc: {}", bits_to_string(up, N),
            // bits_to_string(dnc,N)); lila::Log.out("upc: {} dn : {}",
            // bits_to_string(upc,N), bits_to_string(dn, N));
          }

        for (auto dn : Combinations(N, ndn))
          for (auto upc : Combinations(N - ndn, nup)) {
            auto [up, dnc] = upc_dn_to_up_dnc(upc, dn);
            auto [upc2, dn2] = up_dnc_to_upc_dn(up, dnc);
            REQUIRE(upc == upc2);
            REQUIRE(dn == dn2);
            // lila::Log.out("up : {} dnc: {}", bits_to_string(up, N),
            // bits_to_string(dnc,N)); lila::Log.out("upc: {} dn : {}",
            // bits_to_string(upc,N), bits_to_string(dn, N));
          }

        int n_holes = N - nup - ndn;
        int charge = nup + ndn;
        for (auto holes : Combinations(N, n_holes))
          for (auto spins : Combinations(charge, nup)) {
            auto [up, dn] = holes_spins_to_up_dn(holes, spins, N);
            auto [holes2, spins2] = up_dn_to_holes_spins(up, dn, N);
            // lila::Log.out("holes: {} spins: {}", bits_to_string(holes, N),
            // bits_to_string(spins,charge)); lila::Log.out("up   : {} dn   :
            // {}", bits_to_string(up,N), bits_to_string(dn, N));
            // lila::Log.out("hole2: {} spin2: {}\n", bits_to_string(holes2, N),
            // bits_to_string(spins2,charge));
            REQUIRE(holes2 == holes);
            REQUIRE(spins2 == spins);
          }

        for (auto up : Combinations(N, nup))
          for (auto dn : Combinations(N, ndn)) {
	    if (has_double_occupations(up, dn))
	      continue;
	    
            auto [holes, spins] = up_dn_to_holes_spins(up, dn, N);
            auto [up2, dn2] = holes_spins_to_up_dn(holes, spins, N);
            // lila::Log.out("up   : {} dn   : {}", bits_to_string(up, N),
            //               bits_to_string(dn, N));
            // lila::Log.out("holes: {} spins: {}", bits_to_string(holes, N),
            //               bits_to_string(spins, charge));
            // lila::Log.out("up2: {} dn2: {}\n", bits_to_string(up2, N),
            //               bits_to_string(dn2, charge));
            REQUIRE(up2 == up);
            REQUIRE(dn2 == dn);
          }
      }
}
