#pragma once

#include <xdiag/basis/electron/apply/apply_terms.hpp>
#include <xdiag/common.hpp>

namespace xdiag::basis::electron {

template <typename coeff_t, class Fill>
inline void dispatch(OpSum const &ops, Electron const &block_in,
                     Electron const &block_out, Fill &&fill) try {
  auto const &basis_in = block_in.basis();
  auto const &basis_out = block_out.basis();

  std::visit(
      overload{
          // uint32_t
          [&](BasisNp<uint32_t> const &idx_in,
              BasisNp<uint32_t> const &idx_out) {
            apply_terms<uint32_t, coeff_t, false>(ops, idx_in, idx_out, fill);
          },
          [&](BasisNoNp<uint32_t> const &idx_in,
              BasisNoNp<uint32_t> const &idx_out) {
            apply_terms<uint32_t, coeff_t, false>(ops, idx_in, idx_out, fill);
          },
          [&](BasisSymmetricNp<uint32_t> const &idx_in,
              BasisSymmetricNp<uint32_t> const &idx_out) {
            apply_terms<uint32_t, coeff_t, true>(ops, idx_in, idx_out, fill);
          },
          [&](BasisSymmetricNoNp<uint32_t> const &idx_in,
              BasisSymmetricNoNp<uint32_t> const &idx_out) {
            apply_terms<uint32_t, coeff_t, true>(ops, idx_in, idx_out, fill);
          },

          // uint64_t
          [&](BasisNp<uint64_t> const &idx_in,
              BasisNp<uint64_t> const &idx_out) {
            apply_terms<uint64_t, coeff_t, false>(ops, idx_in, idx_out, fill);
          },
          [&](BasisNoNp<uint64_t> const &idx_in,
              BasisNoNp<uint64_t> const &idx_out) {
            apply_terms<uint64_t, coeff_t, false>(ops, idx_in, idx_out, fill);
          },
          [&](BasisSymmetricNp<uint64_t> const &idx_in,
              BasisSymmetricNp<uint64_t> const &idx_out) {
            apply_terms<uint64_t, coeff_t, true>(ops, idx_in, idx_out, fill);
          },
          [&](BasisSymmetricNoNp<uint64_t> const &idx_in,
              BasisSymmetricNoNp<uint64_t> const &idx_out) {
            apply_terms<uint64_t, coeff_t, true>(ops, idx_in, idx_out, fill);
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

} // namespace xdiag::basis::electron
