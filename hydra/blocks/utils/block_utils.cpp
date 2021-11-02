#include "block_utils.h"

namespace hydra::utils {

bool coupling_is_zero(Bond const &bond, Couplings const &couplings) {
  std::string coupling = bond.coupling();
  return (!couplings.defined(coupling)) ||
         lila::close(couplings[coupling], (complex)0.);
}

bool coupling_is_non_zero(Bond const &bond, Couplings const &couplings) {
  return !coupling_is_zero(bond, couplings);
}

void check_nup_spinhalf(int n_sites, int nup, std::string block_name) {
  if (n_sites < 0)
    lila::Log.err("Error creating {}: n_sites < 0", block_name);
  if ((nup < 0) || (nup > n_sites))
    lila::Log.err("Error creating {}: invalid value of nup", block_name);
}

void check_nup_ndn_electron(int n_sites, int nup, int ndn,
                            std::string block_name) {
  if (n_sites < 0)
    lila::Log.err("Error creating {}: n_sites < 0", block_name);
  if ((nup < 0) || (nup > n_sites))
    lila::Log.err("Error creating {}: invalid value of nup", block_name);
  if ((ndn < 0) || (ndn > n_sites))
    lila::Log.err("Error creating {}: invalid value of ndn", block_name);
}

void check_nup_ndn_tj(int n_sites, int nup, int ndn, std::string block_name) {
  if ((nup < 0) || (nup > n_sites))
    lila::Log.err("Error creating {}: invalid value of nup", block_name);
  if ((ndn < 0) || (ndn > n_sites))
    lila::Log.err("Error creating {}: invalid value of ndn", block_name);
  if (nup + ndn > n_sites)
    lila::Log.err("Error creating {}: invalid value of nup+ndn", block_name);
}

void check_operator_real(BondList const &bonds, Couplings const &cpls,
                         std::string errmsg) {
  if (is_complex(bonds))
    lila::Log.err("Error: cannot {} from complex bonds!", errmsg);
  if (is_complex(cpls))
    lila::Log.err("Error: cannot {} real matrix from complex couplings!",
                  errmsg);
}

void check_symmetric_operator_real(BondList const &bonds, Couplings const &cpls,
                                   Representation const &irrep_in,
                                   Representation const &irrep_out,
                                   std::string errmsg) {
  check_operator_real(bonds, cpls, errmsg);
  if (is_complex(irrep_in) || is_complex(irrep_out))
    lila::Log.err("Error: cannot {} real matrix from complex representation!",
                  errmsg);
}
} // namespace hydra::utils
