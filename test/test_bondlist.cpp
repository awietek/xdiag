#include "catch.hpp"

#include <iostream>

#include "bondlist.h"


TEST_CASE( "bondlist test", "[operators/bondlist]" ) {
  using namespace hydra::operators;

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


  for (auto bond : bl)
    std::cout << bond << std::endl;
  std::cout << std::endl;
  for (auto bond : bl.bonds_of_coupling("J1"))
    std::cout << bond << std::endl;
  std::cout << std::endl;
  for (auto bond : bl.bonds_of_coupling("J2"))
    std::cout << bond << std::endl;

}
