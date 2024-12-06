#pragma once
#include <xdiag/common.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag {

XDIAG_API bool isreal(Op const &op);
XDIAG_API bool isreal(OpSum const &ops);

} // namespace xdiag
