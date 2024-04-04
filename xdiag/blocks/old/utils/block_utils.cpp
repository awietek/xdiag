#include "block_utils.hpp"

#include <xdiag/utils/logger.hpp>
#include <set>
namespace xdiag::utils {

void check_nup_spinhalf(int n_sites, int nup, std::string block_name) {
  if (n_sites < 0)
    Log.err("Error creating {}: n_sites < 0", block_name);
  if ((nup < 0) || (nup > n_sites))
    Log.err("Error creating {}: invalid value of nup", block_name);
}

void check_nup_ndn_electron(int n_sites, int nup, int ndn,
                            std::string block_name) {
  if (n_sites < 0)
    Log.err("Error creating {}: n_sites < 0", block_name);
  if ((nup < 0) || (nup > n_sites))
    Log.err("Error creating {}: invalid value of nup", block_name);
  if ((ndn < 0) || (ndn > n_sites))
    Log.err("Error creating {}: invalid value of ndn", block_name);
}

void check_nup_ndn_tj(int n_sites, int nup, int ndn, std::string block_name) {
  if ((nup < 0) || (nup > n_sites))
    Log.err("Error creating {}: invalid value of nup", block_name);
  if ((ndn < 0) || (ndn > n_sites))
    Log.err("Error creating {}: invalid value of ndn", block_name);
  if (nup + ndn > n_sites)
    Log.err("Error creating {}: invalid value of nup+ndn", block_name);
}

void check_n_sites(int n_sites, PermutationGroup const &permutation_group) {
  if (n_sites != permutation_group.n_sites()) {
    Log.err("Error: number of sites in permutation group is not n_sites!");
  }
}

template <class coeff_t>
std::vector<coeff_t> characters(Representation const &irrep) {
  if constexpr (is_complex<coeff_t>()) {
    return irrep.characters();
  } else {
    return irrep.characters_real();
  }
}

template std::vector<double> characters<double>(Representation const &irrep);
template std::vector<complex> characters<complex>(Representation const &irrep);

} // namespace xdiag::utils
