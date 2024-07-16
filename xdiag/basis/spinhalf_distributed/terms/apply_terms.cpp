#include "apply_terms.hpp"

namespace xdiag::basis::spinhalf_distributed {

template <class basis_t, typename coeff_t>
void apply_terms(OpSum const &ops, basis_t const &basis_in,
                 arma::Col<coeff_t> const &vec_in, basis_t const &basus_out,
                 arma::Col<coeff_t> &vec_out) try {

  for (auto op : ops) {

    if (op.type_defined()) {
      if (op.type() == "EXCHANGE") {
        apply_exchange(op, basis_in, vec_in, basis_out, vec_out);
      } else if (op.type() == "ISING") {
        apply_ising(op, basis_in, vec_in, basis_out, vec_out);
      } else if (op.type() == "SZ") {
        apply_sz(op, basis_in, vec_in, basis_out, vec_out);
      } else if ((op.type() == "S+") || (op.type() == "S-")) {
        apply_spsm(op, basis_in, vec_in, basis_out, vec_out);
      } else if (op.type() == "SCALARCHIRALITY") {
        apply_scalar_chirality(op, basis_in, vec_in, basis_out, vec_out);
      } else {
        XDIAG_THROW(fmt::format("Unknown Op type \"{}\"", op.type()));
      }
    } else {
      apply_non_branching(op, basis_in, vec_in, basis_out, vec_out);
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template void apply_terms(OpSum const &, BasisSz<uint32_t> const &,
                          arma::Col<double> const &, BasisSz<uint32_t> const &,
                          arma::Col<double> &);
template void apply_terms(OpSum const &, BasisSz<uint32_t> const &,
                          arma::Col<complex> const &, BasisSz<uint32_t> const &,
                          arma::Col<complex> &);
template void apply_terms(OpSum const &, BasisSz<uint64_t> const &,
                          arma::Col<double> const &, BasisSz<uint64_t> const &,
                          arma::Col<double> &);
template void apply_terms(OpSum const &, BasisSz<uint64_t> const &,
                          arma::Col<complex> const &, BasisSz<uint64_t> const &,
                          arma::Col<complex> &)

} // namespace xdiag::basis::spinhalf_distributed
