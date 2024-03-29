#pragma once

#include <hydra/basis/basis.h>
#include <hydra/blocks/electron/terms/apply_terms.h>
#include <hydra/common.h>
#include <hydra/utils/logger.h>

namespace hydra::electron {

template <typename coeff_t, class Fill>
inline void dispatch(BondList const &bonds, Electron const &block_in,
                     Electron const &block_out, Fill &&fill) try {
  using namespace basis::electron;

  auto const &basis_in = block_in.basis();
  auto const &basis_out = block_out.basis();

  std::visit(
      overload{
          // uint16_t
          [&](BasisNp<uint16_t> const &idx_in,
              BasisNp<uint16_t> const &idx_out) {
            apply_terms<uint16_t, coeff_t, false>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisNoNp<uint16_t> const &idx_in,
              BasisNoNp<uint16_t> const &idx_out) {
            apply_terms<uint16_t, coeff_t, false>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisSymmetricNp<uint16_t> const &idx_in,
              BasisSymmetricNp<uint16_t> const &idx_out) {
            apply_terms<uint16_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisSymmetricNoNp<uint16_t> const &idx_in,
              BasisSymmetricNoNp<uint16_t> const &idx_out) {
            apply_terms<uint16_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },

          // uint32_t
          [&](BasisNp<uint32_t> const &idx_in,
              BasisNp<uint32_t> const &idx_out) {
            apply_terms<uint32_t, coeff_t, false>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisNoNp<uint32_t> const &idx_in,
              BasisNoNp<uint32_t> const &idx_out) {
            apply_terms<uint32_t, coeff_t, false>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisSymmetricNp<uint32_t> const &idx_in,
              BasisSymmetricNp<uint32_t> const &idx_out) {
            apply_terms<uint32_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisSymmetricNoNp<uint32_t> const &idx_in,
              BasisSymmetricNoNp<uint32_t> const &idx_out) {
            apply_terms<uint32_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },

          // uint64_t
          [&](BasisNp<uint64_t> const &idx_in,
              BasisNp<uint64_t> const &idx_out) {
            apply_terms<uint64_t, coeff_t, false>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisNoNp<uint64_t> const &idx_in,
              BasisNoNp<uint64_t> const &idx_out) {
            apply_terms<uint64_t, coeff_t, false>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisSymmetricNp<uint64_t> const &idx_in,
              BasisSymmetricNp<uint64_t> const &idx_out) {
            apply_terms<uint64_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](BasisSymmetricNoNp<uint64_t> const &idx_in,
              BasisSymmetricNoNp<uint64_t> const &idx_out) {
            apply_terms<uint64_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](auto const &idx_in, auto const &idx_out) {
            HydraThrow(std::logic_error,
                       "Invalid basis or combination of bases");
            (void)idx_in;
            (void)idx_out;
          }},
      basis_in, basis_out);
} catch (...) {
  HydraRethrow("Unable to apply terms on Electron block");
}

} // namespace hydra::electron
