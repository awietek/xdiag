// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "compile.hpp"

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/normal_order.hpp>
#include <xdiag/operators/hc.hpp>
#include <xdiag/operators/valid.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::kernels {

CompiledOpSum compile(OpSum const &ops, Block const &block_in,
                      Block const &block_out) try {
  operators::check_valid(ops);
  check_blocks_match(ops, block_in, block_out);
  if (isdistributed(block_in)) {
    return {ops, block_in, block_out};
  }
  algebra::Algebra alg = algebra::implementation_algebra(block_in);
  return {algebra::normal_order(hc(ops).plain(), alg), block_in, block_out};
}
XDIAG_CATCH

} // namespace xdiag::kernels
