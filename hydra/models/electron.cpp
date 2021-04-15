#include "electron.h"

namespace hydra {
  
  Electron::Electron(int n_sites, QN qn)
    : n_sites_(n_sites), 
      qn_(qn)
  {
    if (qn.defined("Nup") && qn.defined("Ndn"))
      {
	qn_electron qne({qn["Nup"], qn["Ndn"]});
	basis_ = BasisElectron<uint32>(n_sites, qne);
	index_ = IndexElectron<uint32, int64>(basis_);
      }
    else
      {
	printf("Nup/Ndn not defined in QN for Electron\n");
	exit(EXIT_FAILURE);
      }


  }

}
