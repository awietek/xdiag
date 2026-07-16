// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#ifdef XDIAG_DISTRIBUTED

#include <mpi.h>

namespace xdiag::mpi {

template <class coeff_t> MPI_Datatype datatype();

} // namespace xdiag::mpi

#endif
