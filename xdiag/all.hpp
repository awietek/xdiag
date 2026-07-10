/*
Copyright 2026 Alexander Wietek

Licensed under the Apache License, Version 2.0 (the License);
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an AS IS BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#pragma once

#include <xdiag/armadillo.hpp>
#include <xdiag/algebra/ishermitian.hpp>
#include <xdiag/algebra/symmetrize.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/boson.hpp>
#include <xdiag/blocks/electron.hpp>
#include <xdiag/blocks/fermion.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/blocks/tj.hpp>
#include <xdiag/io/file_toml.hpp>
#include <xdiag/io/read.hpp>
#include <xdiag/io/toml/file_toml_handler.hpp>
#include <xdiag/kernels/apply.hpp>
#include <xdiag/kernels/matrix.hpp>
#include <xdiag/kernels/sparse/apply.hpp>
#include <xdiag/kernels/sparse/coo_matrix.hpp>
#include <xdiag/kernels/sparse/csc_matrix.hpp>
#include <xdiag/kernels/sparse/csr_matrix.hpp>
#include <xdiag/kernels/sparse/sparse_matrix_types.hpp>
#include <xdiag/linalg/lanczos/eigs_lanczos.hpp>
#include <xdiag/linalg/lanczos/eigvals_lanczos.hpp>
#include <xdiag/linalg/lobpcg/eigs_lobpcg.hpp>
#include <xdiag/linalg/sparse_diag.hpp>
#include <xdiag/linalg/time_evolution/evolve_lanczos.hpp>
#include <xdiag/linalg/time_evolution/imaginary_time_evolve.hpp>
#include <xdiag/linalg/time_evolution/time_evolve_expokit.hpp>
#include <xdiag/linalg/time_evolution/time_evolve.hpp>
#include <xdiag/operators/hc.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/apply.hpp>
#include <xdiag/states/correlation_matrix.hpp>
#include <xdiag/states/create_state.hpp>
#include <xdiag/states/dot.hpp>
#include <xdiag/states/expect.hpp>
#include <xdiag/states/fill.hpp>
#include <xdiag/states/inner.hpp>
#include <xdiag/states/norm.hpp>
#include <xdiag/states/product_state.hpp>
#include <xdiag/states/random_state.hpp>
#include <xdiag/states/state.hpp>
#include <xdiag/symmetries/cyclic_group.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/permutation.hpp>
#include <xdiag/symmetries/representation.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/say_hello.hpp>
#include <xdiag/utils/timing.hpp>
#include <xdiag/utils/xdiag_api.hpp>
#include <xdiag/utils/xdiag_show.hpp>

#ifdef XDIAG_USE_HDF5
#include <xdiag/io/file_h5.hpp>
#include <xdiag/io/hdf5/file_h5_handler.hpp>
#endif

#ifdef XDIAG_DISTRIBUTED
#include <xdiag/blocks/distributed/electron_distributed.hpp>
#include <xdiag/blocks/distributed/spinhalf_distributed.hpp>
#include <xdiag/blocks/distributed/tj_distributed.hpp>
#endif

#undef XDIAG_API
#undef XDIAG_CATCH
#undef XDIAG_FILL
#undef XDIAG_INSTANTIATE_APPLY
#undef XDIAG_INSTANTIATE_COO_FILL
#undef XDIAG_INSTANTIATE_COO_NNZ
#undef XDIAG_INSTANTIATE_CSR_FILL
#undef XDIAG_INSTANTIATE_CSR_NNZ
#undef XDIAG_INSTANTIATE_KERNELS
#undef XDIAG_INSTANTIATE_KERNELS_BITARRAY
#undef XDIAG_INSTANTIATE_KERNELS_BITARRAY_LONG
#undef XDIAG_INSTANTIATE_MATRIX
#undef XDIAG_LIKELY
#undef XDIAG_OFFSET
#undef XDIAG_RETHROW
#undef XDIAG_SUBLATTICE_UNSTABLE
#undef XDIAG_THROW
#undef XDIAG_TRY_CATCH
#undef XDIAG_UNLIKELY

// exported macros:
// XDIAG_COMPILEDBY
// XDIAG_DIRECTORY
// XDIAG_HOSTNAME
// XDIAG_VERSION
