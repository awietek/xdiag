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
#include <xdiag/misc/api.md>
#include <xdiag/misc/generate_all_hpp.sh>
#include <xdiag/misc/generate_api_md.py>
#include <xdiag/xdiag/algebra/symmetrize.hpp>
#include <xdiag/xdiag/blocks/blocks.hpp>
#include <xdiag/xdiag/blocks/boson.hpp>
#include <xdiag/xdiag/blocks/electron.hpp>
#include <xdiag/xdiag/blocks/fermion.hpp>
#include <xdiag/xdiag/blocks/spinhalf.hpp>
#include <xdiag/xdiag/blocks/tj.hpp>
#include <xdiag/xdiag/io/file_h5.hpp>
#include <xdiag/xdiag/io/file_toml.hpp>
#include <xdiag/xdiag/io/hdf5/file_h5_handler.hpp>
#include <xdiag/xdiag/io/read.hpp>
#include <xdiag/xdiag/kernels/apply.hpp>
#include <xdiag/xdiag/kernels/matrix.hpp>
#include <xdiag/xdiag/kernels/sparse/apply.hpp>
#include <xdiag/xdiag/kernels/sparse/coo_matrix.hpp>
#include <xdiag/xdiag/kernels/sparse/csc_matrix.hpp>
#include <xdiag/xdiag/kernels/sparse/csr_matrix.hpp>
#include <xdiag/xdiag/linalg/lanczos/eigs_lanczos.hpp>
#include <xdiag/xdiag/linalg/lanczos/eigvals_lanczos.hpp>
#include <xdiag/xdiag/linalg/lobpcg/eigs_lobpcg.hpp>
#include <xdiag/xdiag/linalg/sparse_diag.hpp>
#include <xdiag/xdiag/linalg/time_evolution/evolve_lanczos.hpp>
#include <xdiag/xdiag/linalg/time_evolution/imaginary_time_evolve.hpp>
#include <xdiag/xdiag/linalg/time_evolution/time_evolve_expokit.hpp>
#include <xdiag/xdiag/linalg/time_evolution/time_evolve.hpp>
#include <xdiag/xdiag/math/matrix.hpp>
#include <xdiag/xdiag/math/scalar.hpp>
#include <xdiag/xdiag/math/vector.hpp>
#include <xdiag/xdiag/operators/coeff.hpp>
#include <xdiag/xdiag/operators/collect.hpp>
#include <xdiag/xdiag/operators/hc.hpp>
#include <xdiag/xdiag/operators/monomial.hpp>
#include <xdiag/xdiag/operators/op.hpp>
#include <xdiag/xdiag/operators/opsum.hpp>
#include <xdiag/xdiag/states/apply.hpp>
#include <xdiag/xdiag/states/correlation_matrix.hpp>
#include <xdiag/xdiag/states/create_state.hpp>
#include <xdiag/xdiag/states/dot.hpp>
#include <xdiag/xdiag/states/expect.hpp>
#include <xdiag/xdiag/states/fill.hpp>
#include <xdiag/xdiag/states/inner.hpp>
#include <xdiag/xdiag/states/norm.hpp>
#include <xdiag/xdiag/states/product_state.hpp>
#include <xdiag/xdiag/states/random_state.hpp>
#include <xdiag/xdiag/states/state.hpp>
#include <xdiag/xdiag/symmetries/action/site_permutation.hpp>
#include <xdiag/xdiag/symmetries/permutation_group.hpp>
#include <xdiag/xdiag/symmetries/permutation.hpp>
#include <xdiag/xdiag/symmetries/representation.hpp>
#include <xdiag/xdiag/utils/error.hpp>
#include <xdiag/xdiag/utils/logger.hpp>
#include <xdiag/xdiag/utils/say_hello.hpp>
#include <xdiag/xdiag/utils/timing.hpp>
#include <xdiag/xdiag/utils/xdiag_api.hpp>
#include <xdiag/xdiag/utils/xdiag_show.hpp>

#ifdef XDIAG_DISTRIBUTED
#include <xdiag/xdiag/blocks/distributed/electron_distributed.hpp>
#include <xdiag/xdiag/blocks/distributed/spinhalf_distributed.hpp>
#include <xdiag/xdiag/blocks/distributed/tj_distributed.hpp>
#endif

#undef XDIAG_API
#undef XDIAG_CATCH
#undef XDIAG_DEPRECATED
#undef XDIAG_DEPRECATED_EXPORT
#undef XDIAG_DEPRECATED_NO_EXPORT
#undef XDIAG_EXPORT
#undef XDIAG_EXPORT_H
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
#undef XDIAG_NO_DEPRECATED
#undef XDIAG_NO_EXPORT
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
