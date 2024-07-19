#include "apply.hpp"

#include <xdiag/basis/spinhalf_distributed/basis_sz.hpp>

namespace xdiag::basis::spinhalf_distributed {

template <typename coeff_t, class basis_t>
void apply(OpSum const &ops, basis_t const &basis_in,
           arma::Col<coeff_t> const &vec_in, basis_t const &basis_out,
           arma::Col<coeff_t> &vec_out) {
  Log("Arrived!");

  // Determine diagonal and offdiagonal bonds
  auto isdiagonal = [](Op const &op) {
    return (op.type() == "ISING") || (op.type() == "SZ");
  };
  auto isoffdiagonal = [&](Op const &op) { return !isdiagonal(op); };

  OpSum ops_diagonal;
  std::copy_if(ops.begin(), ops.end(), std::back_inserter(ops_diagonal),
               isdiagonal);

  OpSum ops_offdiagonal;
  std::copy_if(ops.begin(), ops.end(), std::back_inserter(ops_offdiagonal),
               isoffdiagonal);

  // Among the offdiagonal bonds, determine which are acting on prefix, postfix
  // or are mixed
  int64_t n_postfix_bits = basis_in.n_postfix_bits();
  auto isprefix = [&](Op const &op) {
    return std::all_of(op.sites().begin(), op.sites().end(),
                       [&](int64_t s) { return s >= n_postfix_bits; });
  };
  auto ispostfix = [&](Op const &op) {
    return std::all_of(op.sites().begin(), op.sites().end(),
                       [&](int64_t s) { return s < n_postfix_bits; });
  };
  auto ismixed = [&](Op const &op) { return !isprefix(op) && !ispostfix(op); };

  OpSum ops_prefix, ops_postfix, ops_mixed;
  std::copy_if(ops_offdiagonal.begin(), ops_offdiagonal.end(),
               std::back_inserter(ops_prefix), isprefix);
  std::copy_if(ops_offdiagonal.begin(), ops_offdiagonal.end(),
               std::back_inserter(ops_postfix), isprefix);
  std::copy_if(ops_offdiagonal.begin(), ops_offdiagonal.end(),
               std::back_inserter(ops_mixed), ismixed);
}

template void apply(OpSum const &, BasisSz<uint32_t> const &,
                    arma::Col<double> const &, BasisSz<uint32_t> const &,
                    arma::Col<double> &);
template void apply(OpSum const &, BasisSz<uint32_t> const &,
                    arma::Col<complex> const &, BasisSz<uint32_t> const &,
                    arma::Col<complex> &);
template void apply(OpSum const &, BasisSz<uint64_t> const &,
                    arma::Col<double> const &, BasisSz<uint64_t> const &,
                    arma::Col<double> &);
template void apply(OpSum const &, BasisSz<uint64_t> const &,
                    arma::Col<complex> const &, BasisSz<uint64_t> const &,
                    arma::Col<complex> &);

} // namespace xdiag::basis::spinhalf_distributed
