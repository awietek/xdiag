#include "compiler.hpp"

#include <xdiag/utils/print_macro.hpp>

namespace xdiag::operators {
void check_bond(Bond const &bond, int64_t n_sites_total, int64_t n_sites_bond,
                bool disjoint, std::string type) {
  check_bond_in_range(bond, n_sites_total);
  check_bond_has_correct_number_of_sites(bond, n_sites_bond);
  if (disjoint) {
    check_bond_has_disjoint_sites(bond);
  }
  if (type == "number") {
    check_bond_coupling_has_type(bond, "double", "complex");
  } else if (type == "matrix") {
    check_bond_coupling_has_type(bond, "arma::mat", "arma::cx_mat");
  } else {
    check_bond_coupling_has_type(bond, type);
  }
}

void check_bond_in_range(Bond const &bond, int64_t n_sites) try {
  for (int64_t s : bond.sites()) {
    if (s >= n_sites) {
      XDIAG_THROW(fmt::format(
          "bond site number {} exceeds range of number of sites = {}", s,
          n_sites));
    } else if (s < 0) {
      XDIAG_THROW(fmt::format("Bond site number {} found to be < 0", s));
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void check_bonds_in_range(BondList const &bonds, int64_t n_sites) try {
  for (Bond bond : bonds) {
    check_bond_in_range(bond, n_sites);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void check_bond_has_correct_number_of_sites(Bond const &bond, int64_t ns) try {
  if (bond.size() != ns) {
    std::string msg = std::string("bond of type \"") + bond.type() +
                      std::string("\": number of sites given is ") +
                      std::to_string(bond.size()) +
                      std::string(" while expected ") + std::to_string(ns);
    XDIAG_THROW(msg);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void check_bond_has_disjoint_sites(Bond const &bond) try {
  if (!bond.sites_disjoint()) {
    XDIAG_THROW(std::string("bond of type \"") + bond.type() +
                std::string("\": sites are not disjoint as expected."));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void check_bond_coupling_has_type(Bond const &bond, std::string type) try {
  std::string bond_type = bond.coupling_type_string();
  if (!(bond_type == type)) {
    XDIAG_THROW(fmt::format("Coupling of bond is expected to be of type \"{}\" "
                            "but received a type \"{}\"",
                            type, bond_type));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void check_bond_coupling_has_type(Bond const &bond, std::string type1,
                                  std::string type2) try {
  std::string bond_type = bond.coupling_type_string();
  if (!(bond_type == type1) && !(bond_type == type2)) {
    XDIAG_THROW(fmt::format("Coupling of bond is expected to be of type either "
                            "\"{}\" or \"{}\" but received a type \"{}\"",
                            type1, type2, bond_type));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

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
    if (bonds.matrix_defined(type)) {
      return true;
    } else {
      return false;
    }
  }
}

complex coupling(Bond const &bond, BondList const &bonds) try {
  if (!coupling_defined(bond, bonds)) {
    XDIAG_THROW(
        fmt::format("coupling \"{}\" is neither defined by Bond nor BondList",
                    bond.coupling_name()));
  }
  if (bond.coupling_defined()) {
    return bond.coupling();
  } else {
    std::string name = bond.coupling_name();
    return bonds.coupling(name);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return 0.;
}

arma::cx_mat matrix(Bond const &bond, BondList const &bonds) try {
  if (!matrix_defined(bond, bonds)) {
    XDIAG_THROW(
        std::string(
            "Error: the following matrix type is neither defined by Bond nor "
            "BondList") +
        bond.type());
  }
  if (bond.matrix_defined()) {
    return bond.matrix();
  } else {
    std::string name = bond.type();
    return bonds.matrix(name);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return arma::cx_mat();
}

BondList compile_explicit_couplings(BondList const &bonds, double precision,
                                    std::string undefined_behavior) try {
  BondList bonds_compiled;

  for (auto it : bonds.matrices()) {
    bonds_compiled.set_matrix(it.first, it.second);
  }
  for (auto it : bonds.couplings()) {
    bonds_compiled.set_coupling(it.first, it.second);
  }

  for (Bond bond : bonds) {

    if (coupling_defined(bond, bonds)) {
      complex cpl = coupling(bond, bonds);

      if (std::abs(cpl) < precision) {
        continue;
      }

      if (bond.type_defined()) {
        bonds_compiled << Bond(bond.type(), cpl, bond.sites());
      } else {
        bonds_compiled << Bond(bond.matrix(), cpl, bond.sites());
      }
    } else {
      if (undefined_behavior == "error") {
        XDIAG_THROW(std::string("undefined coupling in bond: ") +
                    bond.coupling_name());
      } else if (undefined_behavior == "keep") {
        bonds_compiled << bond;
      } else if (undefined_behavior == "ignore") {
      } else {
        XDIAG_THROW("undefined_behaviour must be one of error/keep/ignore");
      }
    }
  }
  return bonds_compiled;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return BondList();
}

BondList compile_explicit_matrices(BondList const &bonds, double precision,
                                   std::string undefined_behavior) try {
  BondList bonds_compiled;
  for (auto it : bonds.matrices()) {
    bonds_compiled.set_matrix(it.first, it.second);
  }
  for (auto it : bonds.couplings()) {
    bonds_compiled.set_coupling(it.first, it.second);
  }

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
        bonds_compiled << Bond(mat, bond.coupling_name(), bond.sites());
      }
    } else {
      if (undefined_behavior == "error") {
        XDIAG_THROW(std::string("undefined matrix in bond: ") +
                    bond.coupling_name());
      } else if (undefined_behavior == "keep") {
        bonds_compiled << bond;
      } else if (undefined_behavior == "ignore") {
      } else {
        XDIAG_THROW(
            "Error: undefined_behaviour must be one of error/keep/ignore");
      }
    }
  }
  return bonds_compiled;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return BondList();
}

BondList compile_explicit(BondList const &bonds, double precision,
                          std::string undefined_behavior) try {
  BondList bonds_compiled;
  bonds_compiled =
      compile_explicit_couplings(bonds, precision, undefined_behavior);
  bonds_compiled =
      compile_explicit_matrices(bonds_compiled, precision, undefined_behavior);
  return bonds_compiled;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return BondList();
}

} // namespace xdiag::operators
