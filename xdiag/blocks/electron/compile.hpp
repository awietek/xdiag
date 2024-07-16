#pragma once

#include <xdiag/common.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::electron {

OpSum compile(OpSum const &ops, int64_t n_sites, double precision);

} // namespace xdiag::electron
