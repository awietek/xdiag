#include "compiler.hpp"

#include <xdiag/utils/print_macro.hpp>

namespace xdiag::operators {



BondList clean_zeros(BondList const &bonds, double precision) try {
  BondList clean_bonds;
  for (auto const &bond : bonds) {

    // if bond is only defined by string, cannot be cleaned
    if (!bond.isexplicit()) {
      XDIAG_THROW("Cannot clean zero bond since is is not explicit. Maybe call "
                  "\"make_explicit\" first");
    }

    Coupling cpl = bond.coupling();
    if (cpl.is<double>()) {
      double cval = cpl.as<double>();
      if (std::abs(cval) > precision) {
        clean_bonds += bond;
      }
    } else if (cpl.is<complex>()) {
      complex cval = cpl.as<complex>();
      if (std::abs(cval) > precision) {
        clean_bonds += bond;
      }
    }
    // coupling is matrix
    else {
      clean_bonds += bond;
    }
  }
  return clean_bonds;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

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
  if (!sites_disjoint(bond)) {
    XDIAG_THROW(std::string("bond of type \"") + bond.type() +
                std::string("\": sites are not disjoint as expected."));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void check_bond_coupling_has_type(Bond const &bond, std::string type) try {
  std::string bond_type = bond.coupling().type();
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
  std::string bond_type = bond.coupling().type();
  if (!(bond_type == type1) && !(bond_type == type2)) {
    XDIAG_THROW(fmt::format("Coupling of bond is expected to be of type either "
                            "\"{}\" or \"{}\" but received a type \"{}\"",
                            type1, type2, bond_type));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::operators
