#pragma once
#ifdef XDIAG_USE_MPI
#include <variant>
#include <xdiag/basis/tj_distributed/basis_np.hpp>
#include <xdiag/common.hpp>

namespace xdiag::basis {
// clang-format off
using BasistJDistributed =
  std::variant<tj_distributed::BasisNp<uint32_t>,
	       tj_distributed::BasisNp<uint64_t>>;
// clang-format on

int64_t dim(BasistJDistributed const &basis);
int64_t size(BasistJDistributed const &basis);
int64_t size_max(BasistJDistributed const &basis);
int64_t size_min(BasistJDistributed const &basis);

template <typename bit_t> bool has_bit_t(BasistJDistributed const &);

} // namespace xdiag::basis

#endif
