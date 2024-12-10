#include "operators.hpp"

#include <xdiag/common.hpp>
#include <xdiag/io/toml/arma_matrix.hpp>
#include <xdiag/io/toml/value.hpp>

namespace xdiag::io {

Scalar scalar(toml::node const &node) try {
  auto real_node = node.value<double>();
  auto cplx_node = node.as_array();
  if (real_node) {
    return Scalar(*real_node);
  } else if (cplx_node) {
    return Scalar(value<complex>(node));
  } else {
    XDIAG_THROW("Invalid TOML format for Scalar. Must be either a real number "
                "or a length 2 array denoting a complex number");
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Coupling coupling(toml::node const &node) try {
  auto real_node = node.value<double>();
  auto cplx_node = node.as_array();
  auto string_node = node.value<std::string>();
  if (real_node || cplx_node) {
    return Coupling(scalar(node));
  } else if (string_node) {
    return Coupling(*string_node);
  } else {
    XDIAG_THROW("Invalid TOML format for Coupling. Must be either a string, a "
                "real number, or a complex number given by a length 2 array.");
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Matrix matrix(toml::node const &node) try {
  // A real matrix should be a MxN array
  // A complex matrix should be a MxNx2 array
  auto array1 = node.as_array();
  if (array1) {
    if (array1->size() > 0) {

      auto array2 = array1[0].as_array();
      if (array2) {
        if (array2->size() > 0) {
          auto array3 = array2[0].as_array();
          auto val = array2[0].value<double>();
          if (array3) {
            if (array3->size() == 2) {
              return Matrix(arma_matrix<complex>(node));
            } else {
              XDIAG_THROW(
                  "TOML node cannot be converted to Matrix: Ill formed three "
                  "dimensional array. If the third dimension is to denote "
                  "complex numbers, it's length needs to be 2");
            }
          } else if (val) {
            return Matrix(arma_matrix<double>(node));
          } else {
            XDIAG_THROW("TOML node cannot be converted to Matrix: Ill formed "
                        "two dimensional array. Entires must either be "
                        "floating point numbers or length two arrays of "
                        "floating point numbers denoting a complex number.");
          }
        } else {
          XDIAG_THROW("TOML node cannot be converted to Matrix: Ill formed two "
                      "dimensional array. Second dimension should not be empty "
                      "arrays.");
        }
      } else {
        XDIAG_THROW("TOML node cannot be converted to Matrix: Ill formed "
                    "array, needs to be two dimensional (MxN) for real matrix "
                    "or three dimensional (MxNx2) for complex matrix.");
      }

    } else { // first array is empty
      return Matrix();
    }
  } else {
    XDIAG_THROW(
        "TOML node cannot be converted to Matrix. Node is not an array.")
  }

} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

static Op op_from_array(toml::array const &array, int64_t start) try {
  int64_t size = array.size();
  if (size < start + 1) {
    XDIAG_THROW("Invalid size of TOML array. Must be at least of size 1");
  }

  auto type_node = array[start].value<std::string>();
  if (!type_node) {
    XDIAG_THROW("Invalid format of TOML array. First entry must be a string");
  }
  std::string type = *type_node;
  if (size == start + 1) { // Operator only defined by type
    return Op(type);
  }

  int64_t sites_end = size;

  // last entry defines a Matrix
  auto matrix_node = array[size - 1].as_array();
  if (matrix_node) {
    --sites_end;
  }

  // The remaining entries are the sites
  std::vector<int64_t> sites(sites_end - start - 1);
  for (int64_t i = start + 1; i < sites_end; ++i) {
    sites[i - start - 1] = value<int64_t>(array[i]);
  }

  if (matrix_node) {
    Matrix mat = matrix(array[size - 1]);
    return Op(type, sites, mat);
  } else {
    return Op(type, sites);
  }

} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Op op(toml::node const &node) try {
  auto array = node.as_array();
  if (array) {
    return op_from_array(*array, 0);
  } else {
    XDIAG_THROW(
        "Error converting TOML to Op: Input TOML node is not an array.");
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

OpSum opsum(toml::node const &node) try {
  auto array = node.as_array();
  if (array) {
    auto ops = OpSum();
    int64_t size = array->size();
    for (int64_t i = 0; i < size; ++i) {
      auto cpl_op_array = array[i].as_array();
      if (!cpl_op_array || (cpl_op_array->size() < 2)) {
        XDIAG_THROW("Cannot convert TOML to OpSum. Entries in OpSum array most "
                    "be arrays which describe both the coupling and the Op");
      }
      Coupling cpl = coupling(cpl_op_array[0]);
      Op op = op_from_array(*cpl_op_array, 1);
      ops += cpl * op;
    }
    return ops;
  } else {
    XDIAG_THROW(
        "Cannot convert TOML to OpSum. Entries must me at least arrays of "
        "length two containing at least a coupling and a type");
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::io
