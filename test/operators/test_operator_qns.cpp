#include "../catch.hpp"

#include <iostream>

#include <hydra/all.h>


TEST_CASE( "operator_qns", "[operators]" ) {
  using namespace hydra;

  for (auto type : {"HB", "HEISENBERG", "EXCHANGE", "ISING", "SZ"})
  {
    Couplings cpls;
    BondList bl;
    bl << Bond(type, "J", {0, 1});
    bl << Bond(type, "J", {1, 2});  
    bl << Bond(type, "J", {2, 3});  
    cpls["J"] = 1.0;
    REQUIRE(utils::spinhalf_nup(bl, cpls) == 0);
    REQUIRE(utils::tj_nup_ndn(bl, cpls) == std::pair<int,int>{0, 0});
    REQUIRE(utils::electron_nup_ndn(bl, cpls) == std::pair<int,int>{0, 0});
  }

  {
    Couplings cpls;
    BondList bl;
    bl << Bond("S+", "hx", 0);
    cpls["hx"] = 1.0;
    REQUIRE(utils::spinhalf_nup(bl, cpls) == 1);
    REQUIRE(utils::tj_nup_ndn(bl, cpls) == std::pair<int,int>{1, -1});
    REQUIRE(utils::electron_nup_ndn(bl, cpls) == std::pair<int,int>{1, -1});
  }

  {
    Couplings cpls;
    BondList bl;
    bl << Bond("S-", "hx", 0);
    cpls["hx"] = 1.0;
    REQUIRE(utils::spinhalf_nup(bl, cpls) == -1);
    REQUIRE(utils::tj_nup_ndn(bl, cpls) == std::pair<int,int>{-1, 1});
    REQUIRE(utils::electron_nup_ndn(bl, cpls) == std::pair<int,int>{-1, 1});
  }

  for (auto type : {"SX", "SY"})
  {
    Couplings cpls;
    BondList bl;
    bl << Bond(type, "hx", 0);
    cpls["hx"] = 1.0;
    REQUIRE(utils::spinhalf_nup(bl, cpls) == undefined_qn);
    REQUIRE(utils::tj_nup_ndn(bl, cpls) == undefined_qns);
    REQUIRE(utils::electron_nup_ndn(bl, cpls) == undefined_qns);
  }

  {
    Couplings cpls;
    BondList bl;
    bl << Bond("Cdagup", "hx", 0);
    cpls["hx"] = 1.0;
    REQUIRE(utils::tj_nup_ndn(bl, cpls) == std::pair<int,int>{1, 0});
    REQUIRE(utils::electron_nup_ndn(bl, cpls) == std::pair<int,int>{1, 0});
  }

  {
    Couplings cpls;
    BondList bl;
    bl << Bond("Cup", "hx", 0);
    cpls["hx"] = 1.0;
    REQUIRE(utils::tj_nup_ndn(bl, cpls) == std::pair<int,int>{-1, 0});
    REQUIRE(utils::electron_nup_ndn(bl, cpls) == std::pair<int,int>{-1, 0});
  }
  {
    Couplings cpls;
    BondList bl;
    bl << Bond("Cdagdn", "hx", 0);
    cpls["hx"] = 1.0;
    REQUIRE(utils::tj_nup_ndn(bl, cpls) == std::pair<int,int>{0, 1});
    REQUIRE(utils::electron_nup_ndn(bl, cpls) == std::pair<int,int>{0, 1});
  }
  {
    Couplings cpls;
    BondList bl;
    bl << Bond("Cdn", "hx", 0);
    cpls["hx"] = 1.0;
    REQUIRE(utils::tj_nup_ndn(bl, cpls) == std::pair<int,int>{0, -1});
    REQUIRE(utils::electron_nup_ndn(bl, cpls) == std::pair<int,int>{0, -1});
  }

}
