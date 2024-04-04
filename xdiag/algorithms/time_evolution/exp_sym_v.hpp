#pragma once

#include <xdiag/algorithms/lanczos/lanczos.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/bondlist.hpp>
#include <xdiag/utils/timing.hpp>

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
