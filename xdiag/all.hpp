#pragma once

#if __cplusplus < 201703L
#error "XDiag requires at least C++17"
#endif

#ifdef XDIAG_USE_MPI
#include <mpi.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#include <xdiag/parallel/omp/omp_utils.hpp>
// #include <xdiag/symmetries/operations/representative_list_omp.hpp>
#endif

#ifdef XDIAG_USE_HDF5
#define ARMA_USE_HDF5
#endif

#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/extern/clara/clara.hpp>
#include <xdiag/extern/fmt/core.hpp>

#include <xdiag/common.hpp>
#include <xdiag/config.hpp>

#include <xdiag/bits/bitops.hpp>
#include <xdiag/utils/close.hpp>
#include <xdiag/utils/iochecks.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/precision.hpp>
#include <xdiag/utils/say_hello.hpp>
#include <xdiag/utils/type_string.hpp>
#include <xdiag/utils/xdiag_show.hpp>

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/apply.hpp>
#include <xdiag/algebra/matrix.hpp>

// Includes for different block types
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/electron.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/blocks/tj.hpp>

#include <xdiag/states/create_state.hpp>
#include <xdiag/states/fill.hpp>
#include <xdiag/states/gpwf.hpp>
#include <xdiag/states/product_state.hpp>
#include <xdiag/states/random_state.hpp>
#include <xdiag/states/state.hpp>

#include <xdiag/random/hash.hpp>
#include <xdiag/random/hash_functions.hpp>
#include <xdiag/random/random_utils.hpp>

#include <xdiag/combinatorics/binomial.hpp>
#include <xdiag/combinatorics/bit_patterns.hpp>
#include <xdiag/combinatorics/combinations.hpp>
#include <xdiag/combinatorics/combinations_index.hpp>
#include <xdiag/combinatorics/combinations_indexing.hpp>
#include <xdiag/combinatorics/fermi_table.hpp>
#include <xdiag/combinatorics/lin_table.hpp>
#include <xdiag/combinatorics/subsets.hpp>
#include <xdiag/combinatorics/subsets_index.hpp>
#include <xdiag/combinatorics/subsets_indexing.hpp>

#include <xdiag/basis/electron/basis_electron.hpp>
#include <xdiag/basis/spinhalf/basis_spinhalf.hpp>
#include <xdiag/basis/tj/basis_tj.hpp>

#include <xdiag/symmetries/continuous_group.hpp>
#include <xdiag/symmetries/generated_group.hpp>
#include <xdiag/symmetries/group_action/group_action.hpp>
#include <xdiag/symmetries/group_action/group_action_lookup.hpp>
#include <xdiag/symmetries/group_action/group_action_sublattice.hpp>
#include <xdiag/symmetries/group_action/sublattice_stability.hpp>
#include <xdiag/symmetries/operations/fermi_sign.hpp>
#include <xdiag/symmetries/operations/group_action_operations.hpp>
#include <xdiag/symmetries/operations/representative_list.hpp>
#include <xdiag/symmetries/operations/symmetry_operations.hpp>
#include <xdiag/symmetries/permutation.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/qn.hpp>
#include <xdiag/symmetries/representation.hpp>

#include <xdiag/operators/non_branching_op.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/operators/scalar.hpp>
#include <xdiag/operators/logic/symmetrize.hpp>
#include <xdiag/operators/logic/real.hpp>
#include <xdiag/operators/logic/hc.hpp>
#include <xdiag/operators/logic/isapprox.hpp>
#include <xdiag/operators/logic/permute.hpp>
#include <xdiag/operators/logic/qns.hpp>
#include <xdiag/operators/logic/valid.hpp>
#include <xdiag/operators/logic/types.hpp>
#include <xdiag/operators/logic/compilation.hpp>

#include <xdiag/algorithms/norm_estimate.hpp>
#include <xdiag/algorithms/sparse_diag.hpp>

#include <xdiag/algorithms/lanczos/eigs_lanczos.hpp>
#include <xdiag/algorithms/lanczos/eigvals_lanczos.hpp>
#include <xdiag/algorithms/lanczos/lanczos.hpp>
#include <xdiag/algorithms/lanczos/lanczos_convergence.hpp>
#include <xdiag/algorithms/lanczos/lanczos_pro.hpp>
#include <xdiag/algorithms/lanczos/tmatrix.hpp>

#include <xdiag/algorithms/arnoldi/arnoldi.hpp>
#include <xdiag/algorithms/arnoldi/arnoldi_to_disk.hpp>
#include <xdiag/algorithms/gram_schmidt/gram_schmidt.hpp>
#include <xdiag/algorithms/gram_schmidt/orthogonalize.hpp>

#include <xdiag/algorithms/time_evolution/exp_sym_v.hpp>
#include <xdiag/algorithms/time_evolution/pade_matrix_exponential.hpp>
#include <xdiag/algorithms/time_evolution/time_evolution.hpp>
#include <xdiag/algorithms/time_evolution/zahexpv.hpp>

#include <xdiag/io/file_h5.hpp>
#include <xdiag/io/file_toml.hpp>

#ifdef XDIAG_USE_MPI
#include <xdiag/parallel/mpi/allreduce.hpp>
#include <xdiag/parallel/mpi/alltoall.hpp>
#include <xdiag/parallel/mpi/buffer.hpp>
#include <xdiag/parallel/mpi/cdot_distributed.hpp>
#include <xdiag/parallel/mpi/communicator.hpp>
#include <xdiag/parallel/mpi/datatype.hpp>
#include <xdiag/parallel/mpi/timing_mpi.hpp>

#include <xdiag/basis/spinhalf_distributed/basis_spinhalf_distributed.hpp>
#include <xdiag/basis/spinhalf_distributed/basis_sz.hpp>

#include <xdiag/basis/tj_distributed/basis_np.hpp>
#include <xdiag/basis/tj_distributed/basis_tj_distributed.hpp>

#include <xdiag/blocks/spinhalf_distributed.hpp>
#include <xdiag/blocks/tj_distributed.hpp>
#endif

#undef XDIAG_THROW
#undef XDIAG_RETHROW
#undef XDIAG_API
