#pragma once

#include <xdiag/basis/basis.hpp>
#include <xdiag/blocks/tj/terms/apply_terms.hpp>
#include <xdiag/common.hpp>
#include <xdiag/utils/logger.hpp>

namespace xdiag::tj {

template <typename coeff_t, class fill_f>
inline void dispatch(BondList const &bonds, tJ const &block_in,
                     tJ const &block_out, fill_f &&fill) try {
  using namespace basis::tj;
  auto const &basis_in = block_in.basis();
  auto const &basis_out = block_out.basis();

  std::visit(
      overload{
          // uint16_t
          [&](BasisNp<uint16_t> const &idx_in,
              BasisNp<uint16_t> const &idx_out) {
            apply_terms<uint16_t, coeff_t, false>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisSymmetricNp<uint16_t> const &idx_in,
              BasisSymmetricNp<uint16_t> const &idx_out) {
            apply_terms<uint16_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          // uint32_t
          [&](BasisNp<uint32_t> const &idx_in,
              BasisNp<uint32_t> const &idx_out) {
            apply_terms<uint32_t, coeff_t, false>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisSymmetricNp<uint32_t> const &idx_in,
              BasisSymmetricNp<uint32_t> const &idx_out) {
            apply_terms<uint32_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          // uint64_t
          [&](BasisNp<uint64_t> const &idx_in,
              BasisNp<uint64_t> const &idx_out) {
            apply_terms<uint64_t, coeff_t, false>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisSymmetricNp<uint64_t> const &idx_in,
              BasisSymmetricNp<uint64_t> const &idx_out) {
            apply_terms<uint64_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](auto const &idx_in, auto const &idx_out) {
            XDiagThrow(std::logic_error,
                       "Invalid basis or combination of bases");
            (void)idx_in;
            (void)idx_out;
          }},
      basis_in, basis_out);
} catch (...) {
  XDiagRethrow("Unable to apply terms on tJ block");
}

} // namespace xdiag::tj
