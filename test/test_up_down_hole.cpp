#include "catch.hpp"

#include <iostream>
#include <hydra/all.h>

template <class state_t, class index_t>
void test_up_down_hole(int n_sites, int n_upspins, int n_downspins){
  using namespace hydra::all;

  using hs_t = Spinhalf<state_t>;
  {
    hs_t hsup(n_sites, n_upspins);
    hs_t hsholes(n_sites - n_upspins, n_sites - n_upspins - n_downspins);
  
    for (state_t ups : hsup)
      for (state_t holes : hsholes)
	{
	  auto downs = up_hole_to_down(ups, holes);
	  auto holes2 = up_down_to_hole(ups, downs);
	  REQUIRE(holes2 == holes);
	}
  }


  {
    hs_t hsdown(n_sites, n_downspins);
    hs_t hsholes(n_sites - n_downspins, n_sites - n_upspins - n_downspins);
  
    for (state_t downs : hsdown)
      for (state_t holes : hsholes)
	{
	  auto ups = down_hole_to_up(downs, holes);
	  auto holes2 = down_up_to_hole(downs, ups);
	  REQUIRE(holes2 == holes);
	}
  }
  

}

TEST_CASE( "updownhole test", "[updownhole]" ) {
  using namespace hydra::all;
  for (int n_sites = 0; n_sites < 8; ++n_sites)
    for (int n_upspins = 0; n_upspins <= n_sites; ++n_upspins)
      for (int n_downspins = 0; n_downspins <= n_sites - n_upspins; ++n_downspins)
      {
	test_up_down_hole<uint32, uint16>(n_sites, n_upspins, n_downspins);
	test_up_down_hole<uint32, uint32>(n_sites, n_upspins, n_downspins);
	test_up_down_hole<uint32, uint64>(n_sites, n_upspins, n_downspins);
	test_up_down_hole<uint64, uint16>(n_sites, n_upspins, n_downspins);
	test_up_down_hole<uint64, uint32>(n_sites, n_upspins, n_downspins);
	test_up_down_hole<uint64, uint64>(n_sites, n_upspins, n_downspins);
      }
}
