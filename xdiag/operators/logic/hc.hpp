#pragma once
#include <xdiag/common.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag {

XDIAG_API Op hc(Op const &op);
XDIAG_API OpSum hc(OpSum const &ops);

} // namespace xdiag
