// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/common.hpp>

#include <xdiag/blocks/blocks.hpp>

#include <xdiag/basis/electron/apply/apply_terms.hpp>
#include <xdiag/basis/spinhalf/apply/apply_terms.hpp>
#include <xdiag/basis/tj/apply/apply_terms.hpp>

namespace xdiag::algebra {

template <typename coeff_t, class fill_f>
inline void apply_dispatch(OpSum const &ops, Spinhalf const &block_in,
                           Spinhalf const &block_out, fill_f fill) try {
  using namespace basis::spinhalf;

  auto const &basis_in = block_in.basis();
  auto const &basis_out = block_out.basis();

  std::visit(overload{// uint32_t
                      [&](BasisSz<uint32_t> const &idx_in,
                          BasisSz<uint32_t> const &idx_out) {
                        apply_terms<coeff_t, false>(ops, idx_in, idx_out, fill);
                      },
                      [&](BasisNoSz<uint32_t> const &idx_in,
                          BasisNoSz<uint32_t> const &idx_out) {
                        apply_terms<coeff_t, false>(ops, idx_in, idx_out, fill);
                      },
                      [&](BasisSymmetricSz<uint32_t> const &idx_in,
                          BasisSymmetricSz<uint32_t> const &idx_out) {
                        apply_terms<coeff_t, true>(ops, idx_in, idx_out, fill);
                      },
                      [&](BasisSymmetricNoSz<uint32_t> const &idx_in,
                          BasisSymmetricNoSz<uint32_t> const &idx_out) {
                        apply_terms<coeff_t, true>(ops, idx_in, idx_out, fill);
                      },

                      // uint64_t
                      [&](BasisSz<uint64_t> const &idx_in,
                          BasisSz<uint64_t> const &idx_out) {
                        apply_terms<coeff_t, false>(ops, idx_in, idx_out, fill);
                      },
                      [&](BasisNoSz<uint64_t> const &idx_in,
                          BasisNoSz<uint64_t> const &idx_out) {
                        apply_terms<coeff_t, false>(ops, idx_in, idx_out, fill);
                      },
                      [&](BasisSymmetricSz<uint64_t> const &idx_in,
                          BasisSymmetricSz<uint64_t> const &idx_out) {
                        apply_terms<coeff_t, true>(ops, idx_in, idx_out, fill);
                      },
                      [&](BasisSymmetricNoSz<uint64_t> const &idx_in,
                          BasisSymmetricNoSz<uint64_t> const &idx_out) {
                        apply_terms<coeff_t, true>(ops, idx_in, idx_out, fill);
                      },
                      [&](BasisSublattice<uint64_t, 1> const &idx_in,
                          BasisSublattice<uint64_t, 1> const &idx_out) {
                        apply_terms<coeff_t, true>(ops, idx_in, idx_out, fill);
                      },
                      [&](BasisSublattice<uint64_t, 2> const &idx_in,
                          BasisSublattice<uint64_t, 2> const &idx_out) {
                        apply_terms<coeff_t, true>(ops, idx_in, idx_out, fill);
                      },
                      [&](BasisSublattice<uint64_t, 3> const &idx_in,
                          BasisSublattice<uint64_t, 3> const &idx_out) {
                        apply_terms<coeff_t, true>(ops, idx_in, idx_out, fill);
                      },
                      [&](BasisSublattice<uint64_t, 4> const &idx_in,
                          BasisSublattice<uint64_t, 4> const &idx_out) {
                        apply_terms<coeff_t, true>(ops, idx_in, idx_out, fill);
                      },
                      [&](BasisSublattice<uint64_t, 5> const &idx_in,
                          BasisSublattice<uint64_t, 5> const &idx_out) {
                        apply_terms<coeff_t, true>(ops, idx_in, idx_out, fill);
                      },

                      [&](auto const &idx_in, auto const &idx_out) {
                        XDIAG_THROW("Invalid basis or combination of bases");
                        (void)idx_in;
                        (void)idx_out;
                      }},
             basis_in, basis_out);
}
XDIAG_CATCH

template <typename coeff_t, class fill_f>
inline void apply_dispatch(OpSum const &ops, tJ const &block_in,
                           tJ const &block_out, fill_f fill) try {
  using namespace basis::tj;

  auto const &basis_in = block_in.basis();
  auto const &basis_out = block_out.basis();

  std::visit(overload{// uint32_t
                      [&](BasisNp<uint32_t> const &idx_in,
                          BasisNp<uint32_t> const &idx_out) {
                        apply_terms<coeff_t, false>(ops, idx_in, idx_out, fill);
                      },
                      [&](BasisSymmetricNp<uint32_t> const &idx_in,
                          BasisSymmetricNp<uint32_t> const &idx_out) {
                        apply_terms<coeff_t, true>(ops, idx_in, idx_out, fill);
                      },
                      // uint64_t
                      [&](BasisNp<uint64_t> const &idx_in,
                          BasisNp<uint64_t> const &idx_out) {
                        apply_terms<coeff_t, false>(ops, idx_in, idx_out, fill);
                      },
                      [&](BasisSymmetricNp<uint64_t> const &idx_in,
                          BasisSymmetricNp<uint64_t> const &idx_out) {
                        apply_terms<coeff_t, true>(ops, idx_in, idx_out, fill);
                      },
                      [&](auto const &idx_in, auto const &idx_out) {
                        XDIAG_THROW("Invalid basis or combination of bases");
                      }},
             basis_in, basis_out);
}
XDIAG_CATCH

template <typename coeff_t, class fill_f>
inline void apply_dispatch(OpSum const &ops, Electron const &block_in,
                           Electron const &block_out, fill_f fill) try {
  using namespace basis::electron;

  auto const &basis_in = block_in.basis();
  auto const &basis_out = block_out.basis();

  std::visit(overload{// uint32_t
                      [&](BasisNp<uint32_t> const &idx_in,
                          BasisNp<uint32_t> const &idx_out) {
                        apply_terms<coeff_t, false>(ops, idx_in, idx_out, fill);
                      },
                      [&](BasisNoNp<uint32_t> const &idx_in,
                          BasisNoNp<uint32_t> const &idx_out) {
                        apply_terms<coeff_t, false>(ops, idx_in, idx_out, fill);
                      },
                      [&](BasisSymmetricNp<uint32_t> const &idx_in,
                          BasisSymmetricNp<uint32_t> const &idx_out) {
                        apply_terms<coeff_t, true>(ops, idx_in, idx_out, fill);
                      },
                      [&](BasisSymmetricNoNp<uint32_t> const &idx_in,
                          BasisSymmetricNoNp<uint32_t> const &idx_out) {
                        apply_terms<coeff_t, true>(ops, idx_in, idx_out, fill);
                      },

                      // uint64_t
                      [&](BasisNp<uint64_t> const &idx_in,
                          BasisNp<uint64_t> const &idx_out) {
                        apply_terms<coeff_t, false>(ops, idx_in, idx_out, fill);
                      },
                      [&](BasisNoNp<uint64_t> const &idx_in,
                          BasisNoNp<uint64_t> const &idx_out) {
                        apply_terms<coeff_t, false>(ops, idx_in, idx_out, fill);
                      },
                      [&](BasisSymmetricNp<uint64_t> const &idx_in,
                          BasisSymmetricNp<uint64_t> const &idx_out) {
                        apply_terms<coeff_t, true>(ops, idx_in, idx_out, fill);
                      },
                      [&](BasisSymmetricNoNp<uint64_t> const &idx_in,
                          BasisSymmetricNoNp<uint64_t> const &idx_out) {
                        apply_terms<coeff_t, true>(ops, idx_in, idx_out, fill);
                      },
                      [&](auto const &idx_in, auto const &idx_out) {
                        XDIAG_THROW("Invalid basis or combination of bases");
                      }},
             basis_in, basis_out);
}
XDIAG_CATCH

#ifdef XDIAG_USE_MPI
template <typename coeff_t>
void apply_dispatch(OpSum const &ops, ElectronDistributed const &block_in,
                    arma::Col<coeff_t> const &vec_in,
                    ElectronDistributed const &block_out,
                    arma::Col<coeff_t> &vec_out);
template <typename coeff_t>
void apply_dispatch(OpSum const &ops, ElectronDistributed const &block_in,
                    arma::Mat<coeff_t> const &vec_in,
                    ElectronDistributed const &block_out,
                    arma::Mat<coeff_t> &vec_out);

template <typename coeff_t>
void apply_dispatch(OpSum const &ops, SpinhalfDistributed const &block_in,
                    arma::Col<coeff_t> const &vec_in,
                    SpinhalfDistributed const &block_out,
                    arma::Col<coeff_t> &vec_out);
template <typename coeff_t>
void apply_dispatch(OpSum const &ops, SpinhalfDistributed const &block_in,
                    arma::Mat<coeff_t> const &vec_in,
                    SpinhalfDistributed const &block_out,
                    arma::Mat<coeff_t> &vec_out);

template <typename coeff_t>
void apply_dispatch(OpSum const &ops, tJDistributed const &block_in,
                    arma::Col<coeff_t> const &vec_in,
                    tJDistributed const &block_out,
                    arma::Col<coeff_t> &vec_out);
template <typename coeff_t>
void apply_dispatch(OpSum const &ops, tJDistributed const &block_in,
                    arma::Mat<coeff_t> const &vec_in,
                    tJDistributed const &block_out,
                    arma::Mat<coeff_t> &vec_out);

#endif

} // namespace xdiag::algebra
