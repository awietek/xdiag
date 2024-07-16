#pragma once
#ifdef XDIAG_USE_MPI

#include <xdiag/basis/basis.hpp>
#include <xdiag/blocks/tj_distributed/terms/apply_terms.hpp>
#include <xdiag/common.hpp>

namespace xdiag::tj_distributed {

template <typename coeff_t>
inline void dispatch(OpSum const &ops, tJDistributed const &block_in,
                     arma::Col<coeff_t> const &vec_in,
                     tJDistributed const &block_out,
                     arma::Col<coeff_t> &vec_out) try {
  using namespace basis::tj_distributed;

  auto const &basis_in = block_in.basis();
  auto const &basis_out = block_out.basis();

  std::visit(overload{// uint16_t
                      [&](BasisNp<uint16_t> const &idx_in,
                          BasisNp<uint16_t> const &idx_out) {
                        apply_terms<uint16_t, coeff_t>(ops, idx_in, vec_in,
                                                       idx_out, vec_out);
                      },
                      [&](BasisNp<uint32_t> const &idx_in,
                          BasisNp<uint32_t> const &idx_out) {
                        apply_terms<uint32_t, coeff_t>(ops, idx_in, vec_in,
                                                       idx_out, vec_out);
                      },
                      [&](BasisNp<uint64_t> const &idx_in,
                          BasisNp<uint64_t> const &idx_out) {
                        apply_terms<uint64_t, coeff_t>(ops, idx_in, vec_in,
                                                       idx_out, vec_out);
                      },
                      [&](auto const &idx_in, auto const &idx_out) {
                        XDIAG_THROW("Invalid basis or combination of bases");
                        (void)idx_in;
                        (void)idx_out;
                      }},
             basis_in, basis_out);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::tj_distributed

#endif
