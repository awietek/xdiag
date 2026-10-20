// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "elements.hpp"

#include <type_traits>

#include <xdiag/kernels/blocks/dispatch_bases.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/variants.hpp>

namespace xdiag::kernels {

// The basis-level kernel: defined in kernels.hpp and explicitly instantiated
// per basis type in kernels/blocks/<block>/kernels.cpp, so this file never
// compiles the term kernels itself.
template <typename block_t, typename coeff_t, typename basis_t>
void elements(OpSum const &ops, basis_t const &rows, basis_t const &cols,
              fill_t<coeff_t> const &fill);

template <typename coeff_t>
void elements(CompiledOpSum const &ops, fill_t<coeff_t> const &fill) try {
  utils::visit_same_type(
      ops.block_out, ops.block_in,
      [&](auto const &bout, auto const &bin) {
        using block_t = std::decay_t<decltype(bout)>;
        if constexpr (is_distributed_v<block_t>) {
          XDIAG_THROW("Cannot build a matrix for a distributed block: its "
                      "Hilbert space is distributed across MPI ranks. Use "
                      "apply(...) instead.");
        } else {
          // ops is the compiled trans(ops): the kernels run it from the rows
          // (block_out) to the columns (block_in), see compile.hpp
          dispatch_basis(bout, bin, [&](auto const &rows, auto const &cols) {
            elements<block_t, coeff_t>(ops.ops, rows, cols, fill);
          });
        }
      },
      "Type mismatch of Block types");
}
XDIAG_CATCH

template void elements(CompiledOpSum const &, fill_t<double> const &);
template void elements(CompiledOpSum const &, fill_t<complex> const &);

} // namespace xdiag::kernels
