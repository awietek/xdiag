#include "non_branching_op.hpp"

#include <tuple>
#include <vector>
#include <xdiag/bits/bitops.hpp>
#include <xdiag/common.hpp>

namespace xdiag::operators {

template <typename coeff_t>
std::vector<arma::Mat<coeff_t>>
decompose_matrix_to_nonbranching(arma::Mat<coeff_t> const &mat,
                                 double precision) try {
  int64_t m = (int64_t)mat.n_rows;
  int64_t n = (int64_t)mat.n_cols;
  if (m != n) {
    XDIAG_THROW("Error: Op matrix is not square");
  }

  std::vector<std::tuple<int64_t, int64_t, coeff_t>> all_entries;

  // Get diagonal elements
  for (int64_t i = 0; i < n; ++i) {
    if (std::abs(mat(i, i)) > precision) {
      all_entries.push_back({i, i, mat(i, i)});
    }
  }

  // Get offidagonal elements
  for (int64_t n_diag = 1; n_diag < n; ++n_diag) {
    for (int64_t i = 0; i < n - n_diag; ++i) {
      if (std::abs(mat(i, i + n_diag)) > precision) {
        all_entries.push_back({i, i + n_diag, mat(i, i + n_diag)});
      }
      if (std::abs(mat(i + n_diag, i)) > precision) {
        all_entries.push_back({i + n_diag, i, mat(i + n_diag, i)});
      }
    }
  }

  // Reduce to minimal number of non-branching terms
  std::vector<arma::Mat<coeff_t>> mats_nb;

  while (all_entries.size() != 0) {
    std::vector<int64_t> forbidden_columns;
    std::vector<int64_t> forbidden_rows;
    std::vector<std::tuple<int64_t, int64_t, coeff_t>> current_entries;
    std::vector<int64_t> delete_entries;
    int64_t i = 0;
    for (auto [row, column, coeff] : all_entries) {

      if ((std::find(forbidden_rows.begin(), forbidden_rows.end(), row) ==
           forbidden_rows.end()) &&
          (std::find(forbidden_columns.begin(), forbidden_columns.end(),
                     column) == forbidden_columns.end())) {
        current_entries.push_back({row, column, coeff});
        forbidden_rows.push_back(row);
        forbidden_columns.push_back(column);
        delete_entries.push_back(i);
      }
      ++i;
    }

    for (int64_t i = delete_entries.size() - 1; i >= 0; --i)
      all_entries.erase(all_entries.begin() + delete_entries[i]);

    // Create non-branching matrix
    arma::Mat<coeff_t> mat_nb(m, n, arma::fill::zeros);
    for (auto [i, j, coeff] : current_entries) {
      mat_nb(i, j) = coeff;
    }
    mats_nb.push_back(mat_nb);
  }
  return mats_nb;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

OpSum non_branching_ops(Op const &op, double precision) try {
  if (op.coupling().is<arma::mat>()) {
    arma::mat mat = op.coupling().as<arma::mat>();
    auto mats_nb = decompose_matrix_to_nonbranching(mat, precision);
    OpSum ops_nb;
    for (auto mat_nb : mats_nb) {
      ops_nb += Op("NONBRANCHINGOP", mat_nb, op.sites());
    }
    return ops_nb;
  } else if (op.coupling().is<arma::cx_mat>()) {
    arma::cx_mat mat = op.coupling().as<arma::cx_mat>();
    auto mats_nb = decompose_matrix_to_nonbranching(mat, precision);
    OpSum ops_nb;
    for (auto mat_nb : mats_nb) {
      ops_nb += Op("NONBRANCHINGOP", mat_nb, op.sites());
    }
    return ops_nb;
  } else {
    XDIAG_THROW(
        "Cannot convert Op to nonbranching OpSum. Coupling is not a matrix");
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

OpSum non_branching_ops(OpSum const &ops, double precision) {
  OpSum ops_nb;
  for (Op const &op : ops) {
    ops_nb = ops_nb + non_branching_ops(op, precision);
  }
  return ops_nb;
}

template <typename coeff_t>
bool is_non_branching_matrix(arma::Mat<coeff_t> const &mat, double precision) {
  for (arma::uword i = 0; i < mat.n_rows; ++i) {
    int64_t non_zero_in_row = 0;
    for (arma::uword j = 0; j < mat.n_cols; ++j) {
      if (std::abs(mat(i, j)) > precision) {
        ++non_zero_in_row;
      }
    }
    if (non_zero_in_row > 1) {
      return false;
    }
  }
  return true;
}

bool is_non_branching_op(Op const &op, double precision) {
  if (op.coupling().is<arma::mat>()) {
    arma::mat mat = op.coupling().as<arma::mat>();
    return is_non_branching_matrix(mat, precision);
  } else if (op.coupling().is<arma::cx_mat>()) {
    arma::cx_mat mat = op.coupling().as<arma::cx_mat>();
    return is_non_branching_matrix(mat, precision);
  } else {
    return false;
  }
}

template <typename bit_t, typename coeff_t>
NonBranchingOp<bit_t, coeff_t>::NonBranchingOp(Op const &op,
                                               double precision) try
    : sites_(op.sites()), dim_(1 << sites_.size()), mask_(0) {
  if (!op.ismatrix()) {
    XDIAG_THROW("Error constructing NonBranchingOp: the coupling of the Op "
                "must be an explicit matrix");
  }

  if (!is_non_branching_op(op, precision)) {
    XDIAG_THROW(
        "trying to create a NonBranchingOp from a Op which is branching");
  }

  for (auto s : op.sites()) {
    mask_ |= ((bit_t)1 << s);
  }
  // int64_t n_sites = 2;
  // Log("X: {}", BSTR(mask_));

  mask_ = ~mask_;
  // Log("Y: {}", BSTR(mask_));

  arma::cx_mat matrix_ = op.coupling().as<arma::cx_mat>();

  // Matrix dimension is 2**(no. sites of op)
  if ((matrix_.n_cols != dim_) || (matrix_.n_rows != dim_)) {
    XDIAG_THROW(fmt::format(
        "Error: invalid matrix dimension for non-branching Op matrix. Expected "
        "dim={}, but received n_rows={} and n_cols={}.",
        dim_, matrix_.n_rows, matrix_.n_cols));
  }

  non_zero_term_ = std::vector<bool>(dim_, false);
  state_applied_ = std::vector<bit_t>(dim_, 0);
  coeff_ = std::vector<coeff_t>(dim_, 0.);

  for (bit_t in = 0; in < dim_; ++in) {
    // int64_t non_zero_in_row = 0;

    for (bit_t out = 0; out < dim_; ++out) {
      if (std::abs(matrix_(out, in)) > precision) {
        non_zero_term_[in] = true;
        state_applied_[in] = out;
        if constexpr (isreal<coeff_t>()) {
          if (std::abs(imag(matrix_(out, in))) > precision) {
            XDIAG_THROW("Error: trying to create a real NonBranchingOp, but "
                        "found a truly complex matrix entry");
          }
          coeff_[in] = real(matrix_(out, in));
        } else {
          coeff_[in] = matrix_(out, in);
        }
        // ++non_zero_in_row;
      }
    }

    // // security check
    // if (non_zero_term_[in]) {
    //   assert(non_zero_in_row == 1);
    // } else {
    //   assert(non_zero_in_row == 0);
    // }
  }

  // int64_t n_sites = 1;
  // for (bit_t in = 0; in < dim_; ++in) {
  //   std::cout << "a " << BSTR(in) << " -> " << BSTR(state_applied_[in]) <<
  //   "
  //   "
  //             << coeff_[in] << " " << non_zero_term_[in] << std::endl;
  // }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename bit_t, typename coeff_t>
bool NonBranchingOp<bit_t, coeff_t>::is_diagonal() const {
  for (bit_t i = 0; i < dim_; ++i) {
    if ((non_zero_term_[i]) && (state_applied_[i] != i)) {
      return false;
    }
  }
  return true;
}

template <typename bit_t, typename coeff_t>
bool NonBranchingOp<bit_t, coeff_t>::non_zero_term(bit_t local_state) const {
  return non_zero_term_[local_state];
}
template <typename bit_t, typename coeff_t>
coeff_t NonBranchingOp<bit_t, coeff_t>::coeff(bit_t local_state) const {
  return coeff_[local_state];
}

template <typename bit_t, typename coeff_t>
std::pair<bit_t, coeff_t>
NonBranchingOp<bit_t, coeff_t>::state_coeff(bit_t local_state) const {
  return {state_applied_[local_state], coeff_[local_state]};
}

template <typename bit_t, typename coeff_t>
bit_t NonBranchingOp<bit_t, coeff_t>::extract_local_state(bit_t state) const {
  bit_t local_state = 0;
  for (int64_t i = 0; i < (int64_t)sites_.size(); ++i) {
    local_state |= bits::gbit(state, sites_[i]) << i;
  }
  return local_state;
}

template <typename bit_t, typename coeff_t>
bit_t NonBranchingOp<bit_t, coeff_t>::deposit_local_state(bit_t local_state,
                                                          bit_t state) const {
  // int64_t n_sites = 2;
  // Log("a: {}", BSTR(state));
  // Log("mask: {}", BSTR(mask_));

  state &= mask_; // clear bits on site
  // Log("b: {}", BSTR(state));
  for (int64_t i = 0; i < (int64_t)sites_.size(); ++i) {
    state |= bits::gbit(local_state, i) << sites_[i];
    // Log("c: {}", BSTR(state));
  }
  return state;
}

template <typename bit_t, typename coeff_t>
int64_t NonBranchingOp<bit_t, coeff_t>::number_difference() const {
  int64_t diff = 0;
  bool first_diff = true;
  for (bit_t state = 0; state < dim_; ++state) {
    if (non_zero_term_[state]) {
      int64_t diff_state =
          bits::popcnt(state_applied_[state]) - bits::popcnt(state);
      if (first_diff) {
        diff = diff_state;
        first_diff = false;
      } else {
        if (diff_state != diff) {
          return undefined;
        }
      }
    }
  }
  return diff;
}

template class NonBranchingOp<uint16_t, double>;
template class NonBranchingOp<uint32_t, double>;
template class NonBranchingOp<uint64_t, double>;
template class NonBranchingOp<uint16_t, complex>;
template class NonBranchingOp<uint32_t, complex>;
template class NonBranchingOp<uint64_t, complex>;

} // namespace xdiag::operators
