#include "type_string.hpp"

#include <xdiag/common.hpp>
#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/symmetries/permutation.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/representation.hpp>

namespace xdiag::utils {
template <> std::string type_string<bool>() { return "bool"; }
template <> std::string type_string<std::string>() { return "std::string"; }

template <> std::string type_string<int8_t>() { return "int8_t"; }
template <> std::string type_string<int16_t>() { return "int16_t"; }
template <> std::string type_string<int32_t>() { return "int32_t"; }
template <> std::string type_string<int64_t>() { return "int64_t"; }

template <> std::string type_string<uint8_t>() { return "uint8_t"; }
template <> std::string type_string<uint16_t>() { return "uint16_t"; }
template <> std::string type_string<uint32_t>() { return "uint32_t"; }
template <> std::string type_string<uint64_t>() { return "uint64_t"; }

template <> std::string type_string<double>() { return "double"; }
template <> std::string type_string<complex>() { return "complex"; }

template <> std::string type_string<std::vector<bool>>() {
  return "std::vector<bool>";
}
template <> std::string type_string<std::vector<std::string>>() {
  return "std::vector<std::string>";
}

template <> std::string type_string<std::vector<int8_t>>() {
  return "std::vector<int8_t>";
}
template <> std::string type_string<std::vector<int16_t>>() {
  return "std::vector<int16_t>";
}
template <> std::string type_string<std::vector<int32_t>>() {
  return "std::vector<int32_t>";
}
template <> std::string type_string<std::vector<int64_t>>() {
  return "std::vector<int64_t>";
}

template <> std::string type_string<std::vector<uint8_t>>() {
  return "std::vector<uint8_t>";
}
template <> std::string type_string<std::vector<uint16_t>>() {
  return "std::vector<uint16_t>";
}
template <> std::string type_string<std::vector<uint32_t>>() {
  return "std::vector<uint32_t>";
}
template <> std::string type_string<std::vector<uint64_t>>() {
  return "std::vector<uint64_t>";
}

template <> std::string type_string<std::vector<double>>() {
  return "std::vector<double>";
}
template <> std::string type_string<std::vector<complex>>() {
  return "std::vector<complex>";
}

template <> std::string type_string<arma::vec>() { return "arma::vec"; }
template <> std::string type_string<arma::cx_vec>() { return "arma::cx_vec"; }
template <> std::string type_string<arma::ivec>() { return "arma::ivec"; }
template <> std::string type_string<arma::uvec>() { return "arma::uvec"; }

template <> std::string type_string<arma::mat>() { return "arma::mat"; }
template <> std::string type_string<arma::cx_mat>() { return "arma::cx_mat"; }
template <> std::string type_string<arma::imat>() { return "arma::imat"; }
template <> std::string type_string<arma::umat>() { return "arma::umat"; }

template <> std::string type_string<Permutation>() { return "Permutation"; }
template <> std::string type_string<PermutationGroup>() {
  return "PermutationGroup";
}
template <> std::string type_string<Representation>() {
  return "Representation";
}
template <> std::string type_string<Op>() { return "Op"; }
template <> std::string type_string<OpSum>() { return "OpSum"; }
} // namespace xdiag::utils
