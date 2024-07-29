#pragma once
#include <variant>

#include <xdiag/basis/tj/basis_np.hpp>
#include <xdiag/basis/tj/basis_symmetric_np.hpp>
#include <xdiag/common.hpp>

namespace xdiag::basis {

// clang-format off
using BasistJ =
  std::variant<tj::BasisNp<uint32_t>,
	       tj::BasisSymmetricNp<uint32_t>,
	       tj::BasisNp<uint64_t>,
	       tj::BasisSymmetricNp<uint64_t>>;
// clang-format on

// clang-format off
using BasistJIterator =
  std::variant<tj::BasisNpIterator<uint32_t>,
	       tj::BasisSymmetricNpIterator<uint32_t>,
	       tj::BasisNpIterator<uint64_t>,
	       tj::BasisSymmetricNpIterator<uint64_t>>;
// clang-format on

int64_t dim(BasistJ const &basis);
int64_t size(BasistJ const &basis);
template <typename bit_t> bool has_bit_t(BasistJ const &);

} // namespace xdiag::basis
