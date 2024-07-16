#pragma once
#include <string>
#include <vector>

#include <xdiag/operators/opsum.hpp>

namespace xdiag::tj {

OpSum compile(OpSum const &ops, int64_t n_sites, double precision);

} // namespace xdiag::tj
