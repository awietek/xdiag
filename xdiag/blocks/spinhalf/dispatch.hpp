#pragma once

#include <xdiag/blocks/spinhalf/terms/apply_terms.hpp>
#include <xdiag/common.hpp>

namespace xdiag::spinhalf {

template <typename coeff_t, class fill_f>
inline void dispatch(OpSum const &ops, Spinhalf const &block_in,
                     Spinhalf const &block_out, fill_f fill,
                     double zero_precision) try {
  using namespace basis::spinhalf;

  auto const &basis_in = block_in.basis();
  auto const &basis_out = block_out.basis();

  std::visit(
      overload{// uint32_t
               [&](BasisSz<uint32_t> const &idx_in,
                   BasisSz<uint32_t> const &idx_out) {
                 apply_terms<uint32_t, coeff_t, false>(ops, idx_in, idx_out,
                                                       fill, zero_precision);
               },
               [&](BasisNoSz<uint32_t> const &idx_in,
                   BasisNoSz<uint32_t> const &idx_out) {
                 apply_terms<uint32_t, coeff_t, false>(ops, idx_in, idx_out,
                                                       fill, zero_precision);
               },
               [&](BasisSymmetricSz<uint32_t> const &idx_in,
                   BasisSymmetricSz<uint32_t> const &idx_out) {
                 apply_terms<uint32_t, coeff_t, true>(ops, idx_in, idx_out,
                                                      fill, zero_precision);
               },
               [&](BasisSymmetricNoSz<uint32_t> const &idx_in,
                   BasisSymmetricNoSz<uint32_t> const &idx_out) {
                 apply_terms<uint32_t, coeff_t, true>(ops, idx_in, idx_out,
                                                      fill, zero_precision);
               },

               // uint64_t
               [&](BasisSz<uint64_t> const &idx_in,
                   BasisSz<uint64_t> const &idx_out) {
                 apply_terms<uint64_t, coeff_t, false>(ops, idx_in, idx_out,
                                                       fill, zero_precision);
               },
               [&](BasisNoSz<uint64_t> const &idx_in,
                   BasisNoSz<uint64_t> const &idx_out) {
                 apply_terms<uint64_t, coeff_t, false>(ops, idx_in, idx_out,
                                                       fill, zero_precision);
               },
               [&](BasisSymmetricSz<uint64_t> const &idx_in,
                   BasisSymmetricSz<uint64_t> const &idx_out) {
                 apply_terms<uint64_t, coeff_t, true>(ops, idx_in, idx_out,
                                                      fill, zero_precision);
               },
               [&](BasisSymmetricNoSz<uint64_t> const &idx_in,
                   BasisSymmetricNoSz<uint64_t> const &idx_out) {
                 apply_terms<uint64_t, coeff_t, true>(ops, idx_in, idx_out,
                                                      fill, zero_precision);
               },
               [&](BasisSublattice<uint64_t, 1> const &idx_in,
                   BasisSublattice<uint64_t, 1> const &idx_out) {
                 apply_terms<uint64_t, coeff_t, true>(ops, idx_in, idx_out,
                                                      fill, zero_precision);
               },
               [&](BasisSublattice<uint64_t, 2> const &idx_in,
                   BasisSublattice<uint64_t, 2> const &idx_out) {
                 apply_terms<uint64_t, coeff_t, true>(ops, idx_in, idx_out,
                                                      fill, zero_precision);
               },
               [&](BasisSublattice<uint64_t, 3> const &idx_in,
                   BasisSublattice<uint64_t, 3> const &idx_out) {
                 apply_terms<uint64_t, coeff_t, true>(ops, idx_in, idx_out,
                                                      fill, zero_precision);
               },
               [&](BasisSublattice<uint64_t, 4> const &idx_in,
                   BasisSublattice<uint64_t, 4> const &idx_out) {
                 apply_terms<uint64_t, coeff_t, true>(ops, idx_in, idx_out,
                                                      fill, zero_precision);
               },
               [&](BasisSublattice<uint64_t, 5> const &idx_in,
                   BasisSublattice<uint64_t, 5> const &idx_out) {
                 apply_terms<uint64_t, coeff_t, true>(ops, idx_in, idx_out,
                                                      fill, zero_precision);
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

} // namespace xdiag::spinhalf
