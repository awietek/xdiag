#include "compile.hpp"
#include <xdiag/operators/compiler.hpp>
#include <xdiag/utils/logger.hpp>

namespace xdiag::electron {

BondList compile(BondList const &bonds, int64_t n_sites) try {
  using namespace operators;
  BondList bonds_explicit = make_explicit(bonds);
  BondList bonds_compiled;
  for (auto bond : bonds_explicit) {
    std::string type = bond.type();

    // Exchange and Ising terms
    if (type == "HB") {
      check_bond(bond, n_sites, 2, true, "number");
      bonds_compiled += Bond("ISING", bond.coupling(), bond.sites());
      bonds_compiled += Bond("EXCHANGE", bond.coupling(), bond.sites());
    } else if (type == "ISING") {
      check_bond(bond, n_sites, 2, true, "number");
      bonds_compiled += Bond("ISING", bond.coupling(), bond.sites());
    } else if (type == "EXCHANGE") {
      check_bond(bond, n_sites, 2, true, "number");
      bonds_compiled += Bond("EXCHANGE", bond.coupling(), bond.sites());

      // Hopping terms
    } else if (type == "HOP") {
      check_bond(bond, n_sites, 2, true, "number");
      bonds_compiled += Bond("HOPUP", bond.coupling(), bond.sites());
      bonds_compiled += Bond("HOPDN", bond.coupling(), bond.sites());
    } else if (type == "HOPUP") {
      check_bond(bond, n_sites, 2, true, "number");
      bonds_compiled += Bond("HOPUP", bond.coupling(), bond.sites());
    } else if (type == "HOPDN") {
      check_bond(bond, n_sites, 2, true, "number");
      bonds_compiled += Bond("HOPDN", bond.coupling(), bond.sites());
    }

    // Number operators
    else if (type == "NUMBER") {
      check_bond(bond, n_sites, 1, false, "number");
      bonds_compiled += Bond("NUMBERUP", bond.coupling(), bond.sites());
      bonds_compiled += Bond("NUMBERDN", bond.coupling(), bond.sites());
    } else if (type == "NUMBERUP") {
      check_bond(bond, n_sites, 1, false, "number");
      bonds_compiled += Bond("NUMBERUP", bond.coupling(), bond.sites());
    } else if (type == "NUMBERDN") {
      check_bond(bond, n_sites, 1, false, "number");
      bonds_compiled += Bond("NUMBERDN", bond.coupling(), bond.sites());
    } else if (type == "SZ") {
      check_bond(bond, n_sites, 1, false, "number");
      if (bond.coupling_is<double>()) {
        bonds_compiled +=
            Bond("NUMBERUP", 0.5 * bond.coupling<double>(), bond.sites());
        bonds_compiled +=
            Bond("NUMBERDN", -0.5 * bond.coupling<double>(), bond.sites());
      } else if (bond.coupling_is<complex>()) {
        bonds_compiled +=
            Bond("NUMBERUP", 0.5 * bond.coupling<complex>(), bond.sites());
        bonds_compiled +=
            Bond("NUMBERDN", -0.5 * bond.coupling<complex>(), bond.sites());
      }
    }

    // Raising Lowering operators
    else if (type == "CDAGUP") {
      check_bond(bond, n_sites, 1, false, "number");
      bonds_compiled += Bond("CDAGUP", bond.coupling(), bond.sites());
    } else if (type == "CDAGDN") {
      check_bond(bond, n_sites, 1, false, "number");
      bonds_compiled += Bond("CDAGDN", bond.coupling(), bond.sites());
    } else if (type == "CUP") {
      check_bond(bond, n_sites, 1, false, "number");
      bonds_compiled += Bond("CUP", bond.coupling(), bond.sites());
    } else if (type == "CDN") {
      check_bond(bond, n_sites, 1, false, "number");
      bonds_compiled += Bond("CDN", bond.coupling(), bond.sites());
    } else {
      XDIAG_THROW(fmt::format("Invalid or undefined Bond type: \"{}\"", type));
    }
  }
  // Set Hubbbbard U term again
  if (bonds.coupling_defined("U")) {
    bonds_compiled["U"] = bonds["U"];
  }

  return bonds_compiled;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return BondList();
}

} // namespace xdiag::electron
