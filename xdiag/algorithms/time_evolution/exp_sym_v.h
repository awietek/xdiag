#pragma once

#include <xdiag/algorithms/lanczos/lanczos.h>
#include <xdiag/common.h>
#include <xdiag/operators/bondlist.h>
#include <xdiag/utils/timing.h>

namespace xdiag {

State exp_sym_v(BondList const &bonds, State state, double tau,
                bool normalize = false, double shift = 0.,
                double precision = 1e-12, int64_t max_iterations = 1000,
                double deflation_tol = 1e-7);

State exp_sym_v(BondList const &bonds, State state, complex tau,
                bool normalize = false, double shift = 0.,
                double precision = 1e-12, int64_t max_iterations = 1000,
                double deflation_tol = 1e-7);

} // namespace xdiag
