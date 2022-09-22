#include "non_branching_bonds.h"

#include <tuple>
#include <vector>

#include <hydra/common.h>

namespace hydra {

BondList non_branching_bonds(Bond const &bond, double precision) {
  if (bond.type_defined()) {
    return BondList({bond});
  }

  arma::cx_mat mat = bond.matrix();
  int m = (int)mat.n_rows;
  int n = (int)mat.n_cols;
  if (m != n) {
    Log.err("Error: bond matrix is not square");
  }

  std::vector<std::tuple<int, int, complex>> all_entries;

  // Get diagonal elements
  for (int i = 0; i < n; ++i) {
    if (std::abs(mat(i, i)) > precision) {
      all_entries.push_back({i, i, mat(i, i)});
    }
  }

  // Get offidagonal elements
  for (int n_diag = 1; n_diag < n; ++n_diag) {
    for (int i = 0; i < n - n_diag; ++i) {
      if (std::abs(mat(i, i + n_diag)) > precision) {
        all_entries.push_back({i, i + n_diag, mat(i, i + n_diag)});
      }
      if (std::abs(mat(i + n_diag, i)) > precision) {
        all_entries.push_back({i + n_diag, i, mat(i + n_diag, i)});
      }
    }
  }

  // Reduce to minimal number of non-branching terms
  std::vector<arma::cx_mat> mats_nb;

  while (all_entries.size() != 0) {
    std::vector<int> forbidden_columns;
    std::vector<int> forbidden_rows;
    std::vector<std::tuple<int, int, complex>> current_entries;
    std::vector<int> delete_entries;
    for (auto [row, column, coeff] : all_entries) {

      if ((std::find(forbidden_rows.begin(), forbidden_rows.end(), row) ==
           forbidden_rows.end()) &&
          (std::find(forbidden_columns.begin(), forbidden_columns.end(),
                     column) == forbidden_columns.end())) {
        current_entries.push_back(all_entries[i]);
        forbidden_rows.push_back(row);
        forbidden_columns.push_back(column);
        delete_entries.push_back(i);
      }
    }

    for (int i = delete_entries.size() - 1; i >= 0; --i)
      all_entries.erase(all_entries.begin() + delete_entries[i]);

    // Create non-branching matrix
    arma::cx_mat mat_nb(m, n, arma::fill::zeros);
    for (auto [i, j, coeff] : current_entries) {
      mat_nb(i, j) = coeff;
    }
    mats_nb.push_back(mat_nb);
  }

  BondList bonds_nb;
  for (arma::cx_mat mat_nb : mats_nb) {
    if (bond.coupling_defined()) {
      bonds_nb << Bond(mat_nb, bond.coupling(), bond.sites());
    } else {
      bonds_nb << Bond(mat_nb, bond.coupling_name(), bond.sites());
    }
  }
  return bonds_nb;
}
BondList non_branching_bonds(BondList const &bonds, double precision) {
  BondList bonds_nb;
  for (Bond bond : bonds) {
    bonds_nb = bonds_nb + non_branching_bonds(bond, precision);
  }
  return bonds_nb;
}

} // namespace hydra
