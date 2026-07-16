// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <cstdint>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::testcases::fermion {

OpSum freefermion_alltoall(int64_t nsites);
OpSum freefermion_alltoall_complex(int64_t nsites);
  
} // namespace xdiag::testcases::fermion
