#include "print.h"

namespace hydra::utils {

  void PrintPretty(const char* identifier, Bond const& bond) 
  {
    printf("%s:\n", identifier);
    printf("%s %s ", bond.type().c_str(), bond.coupling().c_str());
    for (auto site : bond.sites())
      printf("%d ", site);
    printf("\n");      
  }

  void PrintPretty(const char* identifier, BondList const& bondlist) 
  {
    printf("%s:\n", identifier);
    for (auto bond : bondlist)
      {
	printf("%s %s ", bond.type().c_str(), bond.coupling().c_str());
	for (auto site : bond.sites())
	  printf("%d ", site);
	printf("\n");  
      }    
  }

  void PrintPretty(const char* identifier, Couplings const& couplings) 
  {
    printf("%s:\n", identifier);
    for (auto coupling : couplings)
      {
	printf("%s %f %f\n", coupling.first.c_str(), 
	       std::real(coupling.second), std::imag(coupling.second));
      }    
  }

  void PrintPretty(const char* identifier, Tmatrix const& tmat) {
    printf("%s:\n", identifier);
    PrintPretty("alphas", tmat.alphas());
    PrintPretty("betas", tmat.betas());
    PrintPretty("eigenvalues", tmat.eigenvalues());
  }

} // namespace hydra::utils
