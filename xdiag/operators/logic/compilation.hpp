#pragma once

#include <xdiag/operators/opsum.hpp>

namespace xdiag::operators {

OpSum clean_zeros(OpSum const &ops);

OpSum compile_spinhalf(OpSum const &ops);
OpSum compile_tj(OpSum const &ops);
OpSum compile_electron(OpSum const &ops);

template <typename block_t> OpSum compile(OpSum const &ops);

} // namespace xdiag::operators
