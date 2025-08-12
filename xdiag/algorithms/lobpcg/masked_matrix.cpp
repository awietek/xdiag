#include "masked_matrix.hpp"

#include <xdiag/common.hpp>

namespace xdiag::lobpcg {

template <typename coeff_t>
MaskedMatrix<coeff_t>::MaskedMatrix(arma::Mat<coeff_t> const &m,
                                    std::vector<bool> const &mask) try {
  int64_t n_rows = m.n_rows;

  bool current_i_in_range = false;
  int64_t start = 0;
  int64_t end = 0;
  for (int64_t i = 0; i < mask.size(); ++i) {

    if ((!current_i_in_range) && mask[i]) { // found start
      current_i_in_range = true;
      start = i;
    } else if (current_i_in_range && !mask[i]) { // found end -> add matrix
      current_i_in_range = false;
      end = i;
      coeff_t *ptr = m.colptr(start);
      int64_t n_rows = m.n_rows;
      int64_t n_cols = end - start;

      // create lightweigt matrix without copying data
      arma::Mat<coeff_t> m_range = (ptr, n_rows, n_cols, false, true);
      matrices.push_back(m_range);
    }
  }
}
XDIAG_CATCH

template <typename coeff_t> MaskedMatrix<coeff_t>::matrices() const {
  return matrices_;
}

template <typename coeff_t> void MaskedMatrix<coeff_t>::orthonormalize() try {
  for (auto &mat : matrices_) {

  }
 }
XDIAG_CATCH

template <typename coeff_t>
void MaskedMatrix<coeff_t>::orthogonalize_to(arma::Mat<coeff_t> const &X) try {
  for (auto &mat : matrices_) {
    if (X.n_rows != mat.n_rows) {
      XDIAG_THROW(
          fmt::format("Incompatible number of rows. MaskedMatrix: {}, X: {}",
                      mat.n_rows, X.n_rows));
    }
    mat = mat - (X * (X.t() * mat));
  }
}
XDIAG_CATCH

template class MaskedMatrix<double>;
template class MaskedMatrix<complex>;

template <typename coeff_t>
void apply(
    std::function<void(arma::Mat<coeff_t> const &, arma::Mat<coeff_t>)> A,
    MaskedMatrix<coeff_t> const &X, MaskedMatrix<coeff_t> &Y) try {
  auto const &Xmats = X.matrices();
  auto const &Ymats = Y.matrices();
  if (Xmats.size() != Ymats.size()) {
    XDIAG_THROW(
        "The two masked matrices have a different number of submatrices");
  }
  int64_t n_matrices = Xmats.size();
  for (int64_t idx = 0; idx < n_matrices; ++idx) {
    auto const &Xmat = Xmats[idx];
    auto const &Ymat = Ymats[idx];
    if ((Xmat.n_rows != Ymat.n_rows) || (Xmat.n_cols != Ymat.n_cols)) {
      XDIAG_THROW(fmt::format("Two submatrices have incompatible dimensions. "
                              "X: ({},{}), Y: ({},{})",
                              Xmat.n_rows, Xmat.n_cols, Ymat.n_rows,
                              Ymat.n_cols));
    }
    A(Xmat, Ymat);
  }
}
XDIAG_CATCH

template void apply(std::function<void(arma::mat const &, arma::mat)> A,
                    MaskedMatrix<double> const &X, MaskedMatrix<double> &Y);
template void apply(std::function<void(arma::cx_mat const &, arma::cx_mat)> A,
                    MaskedMatrix<complex> const &X, MaskedMatrix<complex> &Y);

} // namespace xdiag::lobpcg
