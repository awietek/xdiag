#include "basis.hpp"

namespace xdiag {

int64_t dim(basis_spinhalf_variant_t const &idxing) {
  return std::visit([&](auto &&idx) { return idx.dim(); }, idxing);
}
int64_t dim(basis_electron_variant_t const &idxing) {
  return std::visit([&](auto &&idx) { return idx.dim(); }, idxing);
}
int64_t dim(basis_tj_variant_t const &idxing) {
  return std::visit([&](auto &&idx) { return idx.dim(); }, idxing);
}

int64_t size(basis_spinhalf_variant_t const &idxing) {
  return std::visit([&](auto &&idx) { return idx.size(); }, idxing);
}
int64_t size(basis_electron_variant_t const &idxing) {
  return std::visit([&](auto &&idx) { return idx.size(); }, idxing);
}
int64_t size(basis_tj_variant_t const &idxing) {
  return std::visit([&](auto &&idx) { return idx.size(); }, idxing);
}

template <typename bit_t>
bool has_bit_t(basis_spinhalf_variant_t const &basis) try {
  using namespace basis::spinhalf;
  using std::is_same;
  return std::visit(
      overload{[](BasisSz<uint16_t> const &) {
                 return is_same<bit_t, uint16_t>::value;
               },
               [](BasisNoSz<uint16_t> const &) {
                 return is_same<bit_t, uint16_t>::value;
               },
               [](BasisSymmetricSz<uint16_t> const &) {
                 return is_same<bit_t, uint16_t>::value;
               },
               [](BasisSymmetricNoSz<uint16_t> const &) {
                 return is_same<bit_t, uint16_t>::value;
               },
               [](BasisSublattice<uint16_t, 1> const &) {
                 return is_same<bit_t, uint16_t>::value;
               },
               [](BasisSublattice<uint16_t, 2> const &) {
                 return is_same<bit_t, uint16_t>::value;
               },
               [](BasisSublattice<uint16_t, 3> const &) {
                 return is_same<bit_t, uint16_t>::value;
               },
               [](BasisSublattice<uint16_t, 4> const &) {
                 return is_same<bit_t, uint16_t>::value;
               },
               [](BasisSublattice<uint16_t, 5> const &) {
                 return is_same<bit_t, uint16_t>::value;
               },
               [](BasisSz<uint32_t> const &) {
                 return is_same<bit_t, uint32_t>::value;
               },
               [](BasisNoSz<uint32_t> const &) {
                 return is_same<bit_t, uint32_t>::value;
               },
               [](BasisSymmetricSz<uint32_t> const &) {
                 return is_same<bit_t, uint32_t>::value;
               },
               [](BasisSymmetricNoSz<uint32_t> const &) {
                 return is_same<bit_t, uint32_t>::value;
               },
               [](BasisSublattice<uint32_t, 1> const &) {
                 return is_same<bit_t, uint32_t>::value;
               },
               [](BasisSublattice<uint32_t, 2> const &) {
                 return is_same<bit_t, uint32_t>::value;
               },
               [](BasisSublattice<uint32_t, 3> const &) {
                 return is_same<bit_t, uint32_t>::value;
               },
               [](BasisSublattice<uint32_t, 4> const &) {
                 return is_same<bit_t, uint32_t>::value;
               },
               [](BasisSublattice<uint32_t, 5> const &) {
                 return is_same<bit_t, uint32_t>::value;
               },
               [](BasisSz<uint64_t> const &) {
                 return is_same<bit_t, uint64_t>::value;
               },
               [](BasisNoSz<uint64_t> const &) {
                 return is_same<bit_t, uint64_t>::value;
               },
               [](BasisSymmetricSz<uint64_t> const &) {
                 return is_same<bit_t, uint64_t>::value;
               },
               [](BasisSymmetricNoSz<uint64_t> const &) {
                 return is_same<bit_t, uint64_t>::value;
               },
               [](BasisSublattice<uint64_t, 1> const &) {
                 return is_same<bit_t, uint64_t>::value;
               },
               [](BasisSublattice<uint64_t, 2> const &) {
                 return is_same<bit_t, uint64_t>::value;
               },
               [](BasisSublattice<uint64_t, 3> const &) {
                 return is_same<bit_t, uint64_t>::value;
               },
               [](BasisSublattice<uint64_t, 4> const &) {
                 return is_same<bit_t, uint64_t>::value;
               },
               [](BasisSublattice<uint64_t, 5> const &) {
                 return is_same<bit_t, uint64_t>::value;
               },
               [](auto &&) { XDIAG_THROW("Invalid Basis type"); }},
      basis);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return false;
}

template bool has_bit_t<uint16_t>(basis_spinhalf_variant_t const &basis);
template bool has_bit_t<uint32_t>(basis_spinhalf_variant_t const &basis);
template bool has_bit_t<uint64_t>(basis_spinhalf_variant_t const &basis);

template <typename bit_t>
bool has_bit_t(basis_electron_variant_t const &basis) try {
  using namespace basis::electron;
  using std::is_same;
  return std::visit(
      overload{[](BasisNp<uint16_t> const &) {
                 return is_same<bit_t, uint16_t>::value;
               },
               [](BasisNoNp<uint16_t> const &) {
                 return is_same<bit_t, uint16_t>::value;
               },
               [](BasisSymmetricNp<uint16_t> const &) {
                 return is_same<bit_t, uint16_t>::value;
               },
               [](BasisSymmetricNoNp<uint16_t> const &) {
                 return is_same<bit_t, uint16_t>::value;
               },
               [](BasisNp<uint32_t> const &) {
                 return is_same<bit_t, uint32_t>::value;
               },
               [](BasisNoNp<uint32_t> const &) {
                 return is_same<bit_t, uint32_t>::value;
               },
               [](BasisSymmetricNp<uint32_t> const &) {
                 return is_same<bit_t, uint32_t>::value;
               },
               [](BasisSymmetricNoNp<uint32_t> const &) {
                 return is_same<bit_t, uint32_t>::value;
               },
               [](BasisNp<uint64_t> const &) {
                 return is_same<bit_t, uint64_t>::value;
               },
               [](BasisNoNp<uint64_t> const &) {
                 return is_same<bit_t, uint64_t>::value;
               },
               [](BasisSymmetricNp<uint64_t> const &) {
                 return is_same<bit_t, uint64_t>::value;
               },
               [](BasisSymmetricNoNp<uint64_t> const &) {
                 return is_same<bit_t, uint64_t>::value;
               },
               [](auto &&) { XDIAG_THROW("Invalid Basis type"); }},
      basis);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return false;
}

template bool has_bit_t<uint16_t>(basis_electron_variant_t const &basis);
template bool has_bit_t<uint32_t>(basis_electron_variant_t const &basis);
template bool has_bit_t<uint64_t>(basis_electron_variant_t const &basis);

template <typename bit_t> bool has_bit_t(basis_tj_variant_t const &basis) try {
  using namespace basis::tj;
  using std::is_same;
  return std::visit(
      overload{[](BasisNp<uint16_t> const &) {
                 return is_same<bit_t, uint16_t>::value;
               },
               [](BasisSymmetricNp<uint16_t> const &) {
                 return is_same<bit_t, uint16_t>::value;
               },
               [](BasisNp<uint32_t> const &) {
                 return is_same<bit_t, uint32_t>::value;
               },
               [](BasisSymmetricNp<uint32_t> const &) {
                 return is_same<bit_t, uint32_t>::value;
               },
               [](BasisNp<uint64_t> const &) {
                 return is_same<bit_t, uint64_t>::value;
               },
               [](BasisSymmetricNp<uint64_t> const &) {
                 return is_same<bit_t, uint64_t>::value;
               },
               [](auto &&) { XDIAG_THROW("Invalid Basis type"); }},
      basis);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return false;
}

template bool has_bit_t<uint16_t>(basis_tj_variant_t const &basis);
template bool has_bit_t<uint32_t>(basis_tj_variant_t const &basis);
template bool has_bit_t<uint64_t>(basis_tj_variant_t const &basis);

#ifdef XDIAG_USE_MPI

int64_t dim(basis_tj_distributed_variant_t const &idxing) {
  return std::visit([](auto &&idx) { return idx.dim(); }, idxing);
}
int64_t size(basis_tj_distributed_variant_t const &idxing) {
  return std::visit([](auto &&idx) { return idx.size(); }, idxing);
}

int64_t size_max(basis_tj_distributed_variant_t const &idxing) {
  return std::visit([](auto &&idx) { return idx.size_max(); }, idxing);
}

int64_t size_min(basis_tj_distributed_variant_t const &idxing) {
  return std::visit([](auto &&idx) { return idx.size_min(); }, idxing);
}

template <typename bit_t>
bool has_bit_t(basis_tj_distributed_variant_t const &basis) try {
  using namespace basis::tj_distributed;
  using std::is_same;
  return std::visit(
      overload{[](BasisNp<uint16_t> const &) {
                 return is_same<bit_t, uint16_t>::value;
               },
               [](BasisNp<uint32_t> const &) {
                 return is_same<bit_t, uint32_t>::value;
               },
               [](BasisNp<uint64_t> const &) {
                 return is_same<bit_t, uint64_t>::value;
               },
               [](auto &&) { XDIAG_THROW("Invalid Basis type"); }},
      basis);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return false;
}

template bool has_bit_t<uint16_t>(basis_tj_distributed_variant_t const &basis);
template bool has_bit_t<uint32_t>(basis_tj_distributed_variant_t const &basis);
template bool has_bit_t<uint64_t>(basis_tj_distributed_variant_t const &basis);

#endif

} // namespace xdiag
