#pragma once
#include <variant>

#include <xdiag/basis/spinhalf/basis_no_sz.h>
#include <xdiag/basis/spinhalf/basis_sublattice.h>
#include <xdiag/basis/spinhalf/basis_symmetric_no_sz.h>
#include <xdiag/basis/spinhalf/basis_symmetric_sz.h>
#include <xdiag/basis/spinhalf/basis_sz.h>

#include <xdiag/basis/tj/basis_np.h>
#include <xdiag/basis/tj/basis_symmetric_np.h>

#include <xdiag/basis/electron/basis_no_np.h>
#include <xdiag/basis/electron/basis_np.h>
#include <xdiag/basis/electron/basis_symmetric_no_np.h>
#include <xdiag/basis/electron/basis_symmetric_np.h>

#ifdef XDIAG_USE_MPI
#include <xdiag/basis/tj_distributed/basis_np.h>
#endif

namespace xdiag {

// clang-format off
using basis_spinhalf_variant_t = std::variant<
    basis::spinhalf::BasisSz<uint16_t>,
    basis::spinhalf::BasisNoSz<uint16_t>,
    basis::spinhalf::BasisSymmetricSz<uint16_t>,
    basis::spinhalf::BasisSymmetricNoSz<uint16_t>,
    basis::spinhalf::BasisSublattice<uint16_t, 1>,
    basis::spinhalf::BasisSublattice<uint16_t, 2>,
    basis::spinhalf::BasisSublattice<uint16_t, 3>,
    basis::spinhalf::BasisSublattice<uint16_t, 4>,
    basis::spinhalf::BasisSublattice<uint16_t, 5>,
    basis::spinhalf::BasisSz<uint32_t>,
    basis::spinhalf::BasisNoSz<uint32_t>,
    basis::spinhalf::BasisSymmetricSz<uint32_t>,
    basis::spinhalf::BasisSymmetricNoSz<uint32_t>,
    basis::spinhalf::BasisSublattice<uint32_t, 1>,
    basis::spinhalf::BasisSublattice<uint32_t, 2>,
    basis::spinhalf::BasisSublattice<uint32_t, 3>,
    basis::spinhalf::BasisSublattice<uint32_t, 4>,
    basis::spinhalf::BasisSublattice<uint32_t, 5>,
    basis::spinhalf::BasisSz<uint64_t>,
    basis::spinhalf::BasisNoSz<uint64_t>,
    basis::spinhalf::BasisSymmetricSz<uint64_t>,
    basis::spinhalf::BasisSymmetricNoSz<uint64_t>,
    basis::spinhalf::BasisSublattice<uint64_t, 1>,
    basis::spinhalf::BasisSublattice<uint64_t, 2>,
    basis::spinhalf::BasisSublattice<uint64_t, 3>,
    basis::spinhalf::BasisSublattice<uint64_t, 4>,
    basis::spinhalf::BasisSublattice<uint64_t, 5>>;
// clang-format on

// clang-format off
using basis_electron_variant_t =
  std::variant<basis::electron::BasisNp<uint16_t>,
	       basis::electron::BasisNoNp<uint16_t>,
	       basis::electron::BasisSymmetricNp<uint16_t>,
	       basis::electron::BasisSymmetricNoNp<uint16_t>,
	       basis::electron::BasisNp<uint32_t>,
	       basis::electron::BasisNoNp<uint32_t>,
	       basis::electron::BasisSymmetricNp<uint32_t>,
	       basis::electron::BasisSymmetricNoNp<uint32_t>,
	       basis::electron::BasisNp<uint64_t>,
	       basis::electron::BasisNoNp<uint64_t>,
	       basis::electron::BasisSymmetricNp<uint64_t>,
	       basis::electron::BasisSymmetricNoNp<uint64_t>>;
// clang-format on

// clang-format off
using basis_tj_variant_t =
  std::variant<basis::tj::BasisNp<uint16_t>,
	       basis::tj::BasisSymmetricNp<uint16_t>,
	       basis::tj::BasisNp<uint32_t>,
	       basis::tj::BasisSymmetricNp<uint32_t>,
	       basis::tj::BasisNp<uint64_t>,
	       basis::tj::BasisSymmetricNp<uint64_t>>;
// clang-format on

int64_t dim(basis_spinhalf_variant_t const &idxing);
int64_t dim(basis_tj_variant_t const &idxing);
int64_t dim(basis_electron_variant_t const &idxing);

int64_t size(basis_spinhalf_variant_t const &idxing);
int64_t size(basis_tj_variant_t const &idxing);
int64_t size(basis_electron_variant_t const &idxing);

template <typename bit_t> bool has_bit_t(basis_spinhalf_variant_t const &);
template <typename bit_t> bool has_bit_t(basis_tj_variant_t const &);
template <typename bit_t> bool has_bit_t(basis_electron_variant_t const &);

template <typename bit_t>
int64_t index(basis_spinhalf_variant_t const &, bit_t spins);
template <typename bit_t>
int64_t index(basis_tj_variant_t const &, bit_t ups, bit_t dns);
template <typename bit_t>
int64_t index(basis_electron_variant_t const &, bit_t ups, bit_t dns);

#ifdef XDIAG_USE_MPI

// clang-format off
using basis_tj_distributed_variant_t =
  std::variant<basis::tj_distributed::BasisNp<uint16_t>,
	       basis::tj_distributed::BasisNp<uint32_t>,
	       basis::tj_distributed::BasisNp<uint64_t>>;
// clang-format on

int64_t dim(basis_tj_distributed_variant_t const &idxing);
int64_t size(basis_tj_distributed_variant_t const &idxing);
int64_t size_max(basis_tj_distributed_variant_t const &idxing);
int64_t size_min(basis_tj_distributed_variant_t const &idxing);

template <typename bit_t>
bool has_bit_t(basis_tj_distributed_variant_t const &);

#endif

} // namespace xdiag
