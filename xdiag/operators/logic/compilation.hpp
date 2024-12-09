#pragma once

#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

#include <xdiag/blocks/electron.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/blocks/tj.hpp>

#ifdef XDIAG_USE_MPI
#include <xdiag/blocks/spinhalf_distributed.hpp>
#include <xdiag/blocks/tj_distributed.hpp>
#endif

namespace xdiag::operators {

OpSum clean_zeros(OpSum const &ops);
OpSum compile_spinhalf(OpSum const &ops);
OpSum compile_tj(OpSum const &ops);
OpSum compile_electron(OpSum const &ops);

} // namespace xdiag::operators
