#include "compiler.h"

namespace hydra::compiler {

bool coupling_defined(Bond const &bond, BondList const &bonds) {
  if (bond.coupling_defined()) {
    return true;
  } else {
    std::string name = bond.coupling_name();
    if (bonds.coupling_defined(name)) {
      return true;
    } else {
      return false;
    }
  }
}

bool matrix_defined(Bond const &bond, BondList const &bonds) {
  if (bond.matrix_defined()) {
    return true;
  } else {
    std::string type = bond.type();
    if (bonds.matrix_defined(name)) {
      return true;
    } else {
      return false;
    }
  }
}

complex coupling(Bond const &bond, BondList const &bonds) {
  if (!coupling_defined(bond, bonds)) {
    Log.err("Error: the following coupling is neither defined by Bond nor "
            "BondList: {}",
            bond.coupling_name());
  }
  if (bond.coupling_defined()) {
    return bond.coupling();
  } else {
    std::string name = bond.coupling_name();
    return bonds.get_coupling(name);
  }
}

arma::cx_mat matrix(Bond const &bond, BondList const &bonds) {
  if (!matrix_defined(bond, bonds)) {
    Log.err("Error: the following matrix type is neither defined by Bond nor "
            "BondList: {}",
            bond.type());
  }
  if (bond.matrix_defined()) {
    return bond.matrix();
  } else {
    std::string name = bond.type();
    return bonds.get_matrix(name);
  }
}

BondList compile_explicit_couplings(BondList const &bonds, double precision,
                                    std::string undefined_behavior) {
  BondList bonds_compiled;
  for (Bond bond : bonds) {

    if (coupling_defined(bond, bonds)) {
      complex cpl = coupling(bond, bonds);

      if (std::abs(cpl) < precision) {
        continue;
      }

      if (bond.type_defined()) {
        bonds_compiled << Bond(bond.type(), cpl, bond.sites());
      } else {
        bonds_compiled << Bond(bond.matrix(), cpl, bond.sites())
      }
    } else {
      if (undefined_behavior == "error") {
        Log.err("Error in compile_explicit_couplings: undefined coupling in "
                "bond: {}",
                bond.coupling_name());
      } else if (undefined_behavior == "keep") {
        bonds_compiled << bond;
      } else if (undefined_behavior == "ignore") {
      } else {
        Log.err("Error in compile_explicit_couplings: undefined_behaviour must "
                "be one of error/keep/ignore");
      }
    }
  }
  return bonds_compiled;
}

BondList compile_explicit_matrices(BondList const &bonds, double precision,
                                   std::string undefined_behavior) {
  BondList bonds_compiled;

  for (Bond bond : bonds) {

    if (matrix_defined(bond, bonds)) {
      arma::cx_mat mat = matrix(bond, bonds);

      // Go through matrix and set small elements to zero
      for (arma::uword j = 0; j < mat.n_cols; ++j) {
        for (arma::uword i = 0; i < mat.n_rows; ++i) {
          if (std::abs(mat(i, j)) < precision) {
            mat(i, j) = 0.;
          }
        }
      }

      if (arma::norm(mat) < precision) {
        continue;
      }

      if (bond.coupling_defined()) {
        bonds_compiled << Bond(mat, bond.coupling(), bond.sites());
      } else {
        bonds_compiled << Bond(mat, bond.coupling_name(), bond.sites())
      }
    } else {
      if (undefined_behavior == "error") {
        Log.err("Error in compile_explicit_matrices: undefined matrix type in "
                "bond: {}",
                bond.type());
      } else if (undefined_behavior == "keep") {
        bonds_compiled << bond;
      } else if (undefined_behavior == "ignore") {
      } else {
        Log.err("Error in compile_explicit_matrices: undefined_behaviour must "
                "be one of error/keep/ignore");
      }
    }
  }
  return bonds_compiled;
}

BondList compile_explicit(BondList const &bonds, double precision,
                          std::string undefined_behavior) {
  BondList bonds_compiled;
  bonds_compiled =
      compile_explicit_couplings(bonds, precision, undefined_behavior);
  bonds_compiled =
      compile_explicit_matrices(bonds_compiled, precision, undefined_behavior);
  return bonds_compiled;
}

} // namespace hydra::compiler
