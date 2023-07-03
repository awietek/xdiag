#include "../catch.hpp"

#include <iostream>

TEST_CASE( "bondlist", "[operators/bondlist]" ) {
  using namespace hydra;

  BondList bl;
  bl << Bond("HB", "J1", {0, 1});
  bl << Bond("HB", "J1", {1, 2});  
  bl << Bond("HB", "J1", {2, 3});  
  bl << Bond("HB", "J1", {3, 4});
  bl << Bond("HB", "J1", {4, 5});   
  bl << Bond("HB", "J1", {5, 0}); 
  bl << Bond("HB", "J2", {0, 2});
  bl << Bond("HB", "J2", {1, 3});  
  bl << Bond("HB", "J2", {2, 4});  
  bl << Bond("HB", "J2", {3, 5});
  bl << Bond("HB", "J2", {4, 0}); 
  bl << Bond("HB", "J2", {5, 1});  
  {
  Bond b1("A", "B", {1, 2, 3});
  Bond b2("A", "B", {2, 3, 4});
  std::vector<int> c = {2, 3};
  REQUIRE(common_sites(b1, b2) == c);
  }
  {
  Bond b1("A", "B", {1, 2, 3});
  Bond b2("A", "B", {6, 7, 8});
  std::vector<int> c = {};
  REQUIRE(common_sites(b1, b2) == c);
  }
  {
  Bond b1("A", "B", {2, 2, 3});
  Bond b2("A", "B", {2, 7, 8});
  std::vector<int> c = {2};
  REQUIRE(common_sites(b1, b2) == c);
  }
  
  // for (auto bond : bl)
  //   std::cout << bond << std::endl;
  // std::cout << std::endl;
  // for (auto bond : bl.bonds_of_coupling("J1"))
  //   std::cout << bond << std::endl;
  // std::cout << std::endl;
  // for (auto bond : bl.bonds_of_coupling("J2"))
  //   std::cout << bond << std::endl;

}
