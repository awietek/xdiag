#include "block_utils.hpp"

#include <xdiag/utils/logger.hpp>
#include <set>
namespace xdiag::utils {

void check_nup_spinhalf(int nsites, int nup, std::string block_name) {
  if (nsites < 0)
    Log.err("Error creating {}: nsites < 0", block_name);
  if ((nup < 0) || (nup > nsites))
    Log.err("Error creating {}: invalid value of nup", block_name);
}

void check_nup_ndn_electron(int nsites, int nup, int ndn,
                            std::string block_name) {
  if (nsites < 0)
    Log.err("Error creating {}: nsites < 0", block_name);
  if ((nup < 0) || (nup > nsites))
    Log.err("Error creating {}: invalid value of nup", block_name);
  if ((ndn < 0) || (ndn > nsites))
    Log.err("Error creating {}: invalid value of ndn", block_name);
}

void check_nup_ndn_tj(int nsites, int nup, int ndn, std::string block_name) {
  if ((nup < 0) || (nup > nsites))
    Log.err("Error creating {}: invalid value of nup", block_name);
  if ((ndn < 0) || (ndn > nsites))
    Log.err("Error creating {}: invalid value of ndn", block_name);
  if (nup + ndn > nsites)
    Log.err("Error creating {}: invalid value of nup+ndn", block_name);
}

void check_nsites(int nsites, PermutationGroup const &permutation_group) {
  if (nsites != permutation_group.nsites()) {
    Log.err("Error: number of sites in permutation group is not nsites!");
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
