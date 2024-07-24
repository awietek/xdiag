#include "apply_terms.hpp"

#include <xdiag/basis/spinhalf_distributed/apply/apply_exchange.hpp>
#include <xdiag/basis/spinhalf_distributed/apply/apply_ising.hpp>
#include <xdiag/basis/spinhalf_distributed/apply/apply_spsm.hpp>
#include <xdiag/basis/spinhalf_distributed/apply/apply_sz.hpp>

#include <xdiag/basis/spinhalf_distributed/basis_sz.hpp>
#include <xdiag/basis/spinhalf_distributed/transpose.hpp>

namespace xdiag::basis::spinhalf_distributed {

template <typename coeff_t, class basis_t>
void apply_terms(OpSum const &ops, basis_t const &basis_in,
                 arma::Col<coeff_t> const &vec_in, basis_t const &basis_out,
                 arma::Col<coeff_t> &vec_out) try {

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
               std::back_inserter(ops_postfix), ispostfix);
  std::copy_if(ops_offdiagonal.begin(), ops_offdiagonal.end(),
               std::back_inserter(ops_mixed), ismixed);

  // Diagonal operators
  for (Op op : ops_diagonal) {
    std::string type = op.type();
    if (type == "ISING") {
      apply_ising(op, basis_in, vec_in, vec_out);
    } else if (type == "SZ") {
      apply_sz(op, basis_in, vec_in, vec_out);
    } else {
      XDIAG_THROW(fmt::format("Unknown bond of type \"{}\"", type));
    }
  }

  //////////////////////////
  // Apply postfix operators
  for (Op op : ops_postfix) {
    std::string type = op.type();
    if (type == "EXCHANGE") {
      apply_exchange_postfix(op, basis_in, vec_in, vec_out);
    } else if ((type == "S+") || (type == "S-")) {
      apply_spsm_postfix(op, basis_in, vec_in, basis_out, vec_out);
    } else {
      XDIAG_THROW(fmt::format("Unknown bond of type \"{}\"", type));
    }
  }

  ////////////////////////////////////////////////////////////////////////////
  // Apply prefix operators (result is computed frmo send_buffer and stored in
  // recv_buffer)
  int64_t buffer_size = std::max(basis_out.size_max(), basis_in.size_max());
  mpi::buffer.reserve<coeff_t>(buffer_size);
  coeff_t *send_buffer = mpi::buffer.send<coeff_t>();
  coeff_t *recv_buffer = mpi::buffer.recv<coeff_t>();

  // Transpose to postfix | prefix order
  transpose(basis_in, vec_in.memptr(), false);

  for (Op op : ops_prefix) {
    std::string type = op.type();
    if (type == "EXCHANGE") {
      apply_exchange_prefix<basis_t, coeff_t>(op, basis_in);
    } else if ((type == "S+") || (type == "S-")) {
      apply_spsm_prefix<basis_t, coeff_t>(op, basis_in, basis_out);
    } else {
      XDIAG_THROW(fmt::format("Unknown bond of type \"{}\"", type));
    }
  }

  // Transpose back to prefix | postfix order
  transpose(basis_out, recv_buffer, true);

  // Fill contents of send buffer into vec_out
  for (int64_t idx = 0; idx < vec_out.size(); ++idx) {
    vec_out(idx) += send_buffer[idx];
  }

  /////////////////////////////
  // apply mixed operators
  for (Op op : ops_mixed) {
    std::string type = op.type();
    if (type == "EXCHANGE") {
      apply_exchange_mixed(op, basis_in, vec_in, vec_out);
    } else {
      XDIAG_THROW(fmt::format("Unknown bond of type \"{}\"", type));
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
                          arma::Col<complex> &);

} // namespace xdiag::basis::spinhalf_distributed
