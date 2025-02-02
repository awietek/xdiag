#pragma once

#include <optional>
#include <string>
#include <vector>

#include <xdiag/common.hpp>
#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/utils/matrix.hpp>

namespace xdiag {

class Op {
public:
  XDIAG_API Op() = default;
  explicit XDIAG_API Op(std::string type);

  XDIAG_API Op(std::string type, int64_t site);
  XDIAG_API Op(std::string type, std::vector<int64_t> const &sites);

  XDIAG_API Op(std::string type, int64_t site, arma::mat const &matrix);
  XDIAG_API Op(std::string type, std::vector<int64_t> const &sites,
               arma::mat const &matrix);

  XDIAG_API Op(std::string type, int64_t site, arma::cx_mat const &matrix);
  XDIAG_API Op(std::string type, std::vector<int64_t> const &sites,
               arma::cx_mat const &matrix);

  Op(std::string type, int64_t site, Matrix const &matrix);
  Op(std::string type, std::vector<int64_t> const &sites, Matrix const &matrix);

  XDIAG_API std::string type() const;
  XDIAG_API int64_t size() const;
  XDIAG_API int64_t operator[](int64_t idx) const;
  XDIAG_API std::vector<int64_t> const &sites() const;

  bool hassites() const;
  bool hasmatrix() const;
  Matrix const &matrix() const;

  XDIAG_API bool operator==(const Op &rhs) const;
  XDIAG_API bool operator!=(const Op &rhs) const;

private:
  std::string type_;
  std::optional<std::vector<int64_t>> sites_;
  std::optional<Matrix> matrix_;
};

XDIAG_API std::ostream &operator<<(std::ostream &out, Op const &op);
XDIAG_API std::string to_string(Op const &op);

// helper for julia wrapper
class VectorOp {
public:
  VectorOp() = default;
  void push_back(Op const &op);
  std::vector<Op> vector() const;

private:
  std::vector<Op> v_;
};

} // namespace xdiag
