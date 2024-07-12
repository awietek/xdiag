#include "apply_terms.hpp"

namespace xdiag::basis::spinhalf_distributed {

template <class basis_t, typename coeff_t>
void apply_terms(BondList const &bonds, basis_t const &basis_in,
                 arma::Col<coeff_t> const &vec_in, basis_t const &basus_out,
                 arma::Col<coeff_t> &vec_out) try {

  for (auto bond : bonds) {

    if (bond.type_defined()) {
      if (bond.type() == "EXCHANGE") {
        apply_exchange(bond, basis_in, vec_in, basis_out, vec_out);
      } else if (bond.type() == "ISING") {
        apply_ising(bond, basis_in, vec_in, basis_out, vec_out);
      } else if (bond.type() == "SZ") {
        apply_sz(bond, basis_in, vec_in, basis_out, vec_out);
      } else if ((bond.type() == "S+") || (bond.type() == "S-")) {
        apply_spsm(bond, basis_in, vec_in, basis_out, vec_out);
      } else if (bond.type() == "SCALARCHIRALITY") {
        apply_scalar_chirality(bond, basis_in, vec_in, basis_out, vec_out);
      } else {
        XDIAG_THROW(fmt::format("Unknown bond type \"{}\"", bond.type()));
      }
    } else {
      apply_non_branching(bond, basis_in, vec_in, basis_out, vec_out);
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template void apply_terms(BondList const &bonds,
                          BasisSz<uint32_t> const &block_in,
                          arma::Col<double> const &vec_in,
                          BasisSz<uint32_t> const &block_out,
                          arma::Col<double> &vec_out);
template void apply_terms(BondList const &bonds,
                          BasisSz<uint32_t> const &block_in,
                          arma::Col<complex> const &vec_in,
                          BasisSz<uint32_t> const &block_out,
                          arma::Col<complex> &vec_out);
template void apply_terms(BondList const &bonds,
                          BasisSz<uint64_t> const &block_in,
                          arma::Col<double> const &vec_in,
                          BasisSz<uint64_t> const &block_out,
                          arma::Col<double> &vec_out);
template void apply_terms(BondList const &bonds,
                          BasisSz<uint64_t> const &block_in,
                          arma::Col<complex> const &vec_in,
                          BasisSz<uint64_t> const &block_out,
                          arma::Col<complex> &vec_out)

} // namespace xdiag::basis::spinhalf_distributed
