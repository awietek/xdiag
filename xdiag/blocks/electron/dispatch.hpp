#pragma once

#include <xdiag/basis/basis.hpp>
#include <xdiag/blocks/electron/terms/apply_terms.hpp>
#include <xdiag/common.hpp>
#include <xdiag/utils/logger.hpp>

namespace xdiag::electron {

template <typename coeff_t, class Fill>
inline void dispatch(OpSum const &ops, Electron const &block_in,
                     Electron const &block_out, Fill &&fill,
                     double zero_precision) try {
  using namespace basis::electron;

  auto const &basis_in = block_in.basis();
  auto const &basis_out = block_out.basis();

  std::visit(
      overload{// uint16_t
               [&](BasisNp<uint16_t> const &idx_in,
                   BasisNp<uint16_t> const &idx_out) {
                 apply_terms<uint16_t, coeff_t, false>(ops, idx_in, idx_out,
                                                       fill, zero_precision);
               },
               [&](BasisNoNp<uint16_t> const &idx_in,
                   BasisNoNp<uint16_t> const &idx_out) {
                 apply_terms<uint16_t, coeff_t, false>(ops, idx_in, idx_out,
                                                       fill, zero_precision);
               },
               [&](BasisSymmetricNp<uint16_t> const &idx_in,
                   BasisSymmetricNp<uint16_t> const &idx_out) {
                 apply_terms<uint16_t, coeff_t, true>(ops, idx_in, idx_out,
                                                      fill, zero_precision);
               },
               [&](BasisSymmetricNoNp<uint16_t> const &idx_in,
                   BasisSymmetricNoNp<uint16_t> const &idx_out) {
                 apply_terms<uint16_t, coeff_t, true>(ops, idx_in, idx_out,
                                                      fill, zero_precision);
               },

               // uint32_t
               [&](BasisNp<uint32_t> const &idx_in,
                   BasisNp<uint32_t> const &idx_out) {
                 apply_terms<uint32_t, coeff_t, false>(ops, idx_in, idx_out,
                                                       fill, zero_precision);
               },
               [&](BasisNoNp<uint32_t> const &idx_in,
                   BasisNoNp<uint32_t> const &idx_out) {
                 apply_terms<uint32_t, coeff_t, false>(ops, idx_in, idx_out,
                                                       fill, zero_precision);
               },
               [&](BasisSymmetricNp<uint32_t> const &idx_in,
                   BasisSymmetricNp<uint32_t> const &idx_out) {
                 apply_terms<uint32_t, coeff_t, true>(ops, idx_in, idx_out,
                                                      fill, zero_precision);
               },
               [&](BasisSymmetricNoNp<uint32_t> const &idx_in,
                   BasisSymmetricNoNp<uint32_t> const &idx_out) {
                 apply_terms<uint32_t, coeff_t, true>(ops, idx_in, idx_out,
                                                      fill, zero_precision);
               },

               // uint64_t
               [&](BasisNp<uint64_t> const &idx_in,
                   BasisNp<uint64_t> const &idx_out) {
                 apply_terms<uint64_t, coeff_t, false>(ops, idx_in, idx_out,
                                                       fill, zero_precision);
               },
               [&](BasisNoNp<uint64_t> const &idx_in,
                   BasisNoNp<uint64_t> const &idx_out) {
                 apply_terms<uint64_t, coeff_t, false>(ops, idx_in, idx_out,
                                                       fill, zero_precision);
               },
               [&](BasisSymmetricNp<uint64_t> const &idx_in,
                   BasisSymmetricNp<uint64_t> const &idx_out) {
                 apply_terms<uint64_t, coeff_t, true>(ops, idx_in, idx_out,
                                                      fill, zero_precision);
               },
               [&](BasisSymmetricNoNp<uint64_t> const &idx_in,
                   BasisSymmetricNoNp<uint64_t> const &idx_out) {
                 apply_terms<uint64_t, coeff_t, true>(ops, idx_in, idx_out,
                                                      fill, zero_precision);
               },
               [&](auto const &idx_in, auto const &idx_out) {
                 XDIAG_THROW("Invalid basis or combination of bases");
                 (void)idx_in;
                 (void)idx_out;
               }},
      basis_in, basis_out);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

} // namespace xdiag::electron
