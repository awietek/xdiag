#include "order.hpp"

#include <algorithm>
#include <xdiag/operators/logic/types.hpp>

namespace xdiag {

std::pair<Scalar, Op> order(Scalar const &alpha, Op const &op);

bool less(Op const &o1, Op const &o2);
OpSum order(OpSum const &ops);

} // namespace xdiag
