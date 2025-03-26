#pragma once

#include <xdiag/common.hpp>

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/apply.hpp>
#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/algorithms/lanczos/eigs_lanczos.hpp>
#include <xdiag/algorithms/lanczos/eigvals_lanczos.hpp>
#include <xdiag/algorithms/sparse_diag.hpp>
#include <xdiag/algorithms/time_evolution/evolve_lanczos.hpp>
#include <xdiag/algorithms/time_evolution/imaginary_time_evolve.hpp>
#include <xdiag/algorithms/time_evolution/time_evolve.hpp>
#include <xdiag/algorithms/time_evolution/time_evolve_expokit.hpp>
#include <xdiag/blocks/electron.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/blocks/tj.hpp>
#include <xdiag/io/file_h5.hpp>
#include <xdiag/io/file_toml.hpp>
#include <xdiag/io/read.hpp>
#include <xdiag/operators/coupling.hpp>
#include <xdiag/operators/logic/block.hpp>
#include <xdiag/operators/logic/hc.hpp>
#include <xdiag/operators/logic/isapprox.hpp>
#include <xdiag/operators/logic/qns.hpp>
#include <xdiag/operators/logic/real.hpp>
#include <xdiag/operators/logic/symmetrize.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/create_state.hpp>
#include <xdiag/states/fill.hpp>
#include <xdiag/states/product_state.hpp>
#include <xdiag/states/random_state.hpp>
#include <xdiag/states/state.hpp>
#include <xdiag/symmetries/permutation.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/representation.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/say_hello.hpp>
#include <xdiag/utils/scalar.hpp>
#include <xdiag/utils/xdiag_api.hpp>
#include <xdiag/utils/xdiag_show.hpp>

#ifdef XDIAG_USE_MPI
#include <xdiag/blocks/spinhalf_distributed.hpp>
#include <xdiag/blocks/tj_distributed.hpp>
#include <xdiag/blocks/electron_distributed.hpp>
#endif

#undef XDIAG_THROW
#undef XDIAG_RETHROW
#undef XDIAG_API
#undef XDIAG_OFFSET
#undef XDIAG_PI
