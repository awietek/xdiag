#pragma once
#ifdef XDIAG_USE_MPI
#include <variant>
#include <xdiag/common.hpp>
#include <xdiag/basis/spinhalf_distributed/basis_sz.hpp>

namespace xdiag::basis {
// clang-format off
using BasisSpinhalfDistributed =
  std::variant<basis::spinhalf_distributed::BasisSz<uint32_t>,
	       basis::spinhalf_distributed::BasisSz<uint64_t>>;
// clang-format on

int64_t dim(BasisSpinhalfDistributed const &basis);
int64_t size(BasisSpinhalfDistributed const &basis);
int64_t size_max(BasisSpinhalfDistributed const &basis);
int64_t size_min(BasisSpinhalfDistributed const &basis);

template <typename bit_t> bool has_bit_t(BasisSpinhalfDistributed const &);

} // namespace xdiag::basis

#endif
