#include "compile.hpp"

#include <xdiag/operators/compiler.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/print_macro.hpp>

namespace xdiag::tj {

BondList compile(BondList const &bonds, int64_t n_sites, double precision) try {
  using namespace operators;
  BondList bonds_explicit = make_explicit(bonds);
  BondList bonds_clean = clean_zeros(bonds_explicit, precision);

  BondList bonds_compiled;
  for (auto bond : bonds_clean) {
    std::string type = bond.type();

    // Exchange and Ising terms
    if (type == "HB") {
      check_bond(bond, n_sites, 2, true, "number");
      bonds_compiled += Bond("ISING", bond.coupling(), bond.sites());
      bonds_compiled += Bond("EXCHANGE", bond.coupling(), bond.sites());
    } else if (type == "TJHB") {
      check_bond(bond, n_sites, 2, true, "number");
      bonds_compiled += Bond("TJISING", bond.coupling(), bond.sites());
      bonds_compiled += Bond("EXCHANGE", bond.coupling(), bond.sites());
    } else if (type == "ISING") {
      check_bond(bond, n_sites, 2, true, "number");
      bonds_compiled += bond;
    } else if (type == "TJISING") {
      check_bond(bond, n_sites, 2, true, "number");
      bonds_compiled += bond;
    } else if (type == "EXCHANGE") {
      check_bond(bond, n_sites, 2, true, "number");
      bonds_compiled += bond;
    }

    // Hopping terms
    else if (type == "HOP") {
      check_bond(bond, n_sites, 2, true, "number");
      bonds_compiled += Bond("HOPUP", bond.coupling(), bond.sites());
      bonds_compiled += Bond("HOPDN", bond.coupling(), bond.sites());
    } else if (type == "HOPUP") {
      check_bond(bond, n_sites, 2, true, "number");
      bonds_compiled += bond;
    } else if (type == "HOPDN") {
      check_bond(bond, n_sites, 2, true, "number");
      bonds_compiled += bond;
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
      if (bond.coupling().is<double>()) {
        bonds_compiled +=
            Bond("NUMBERUP", 0.5 * bond.coupling().as<double>(), bond.sites());
        bonds_compiled +=
            Bond("NUMBERDN", -0.5 * bond.coupling().as<double>(), bond.sites());
      } else if (bond.coupling().is<complex>()) {
        bonds_compiled +=
            Bond("NUMBERUP", 0.5 * bond.coupling().as<complex>(), bond.sites());
        bonds_compiled += Bond("NUMBERDN", -0.5 * bond.coupling().as<complex>(),
                               bond.sites());
      }
    }

    // Raising Lowering operators
    else if (type == "CDAGUP") {
      check_bond(bond, n_sites, 1, false, "number");
      bonds_compiled += bond;
    } else if (type == "CDAGDN") {
      check_bond(bond, n_sites, 1, false, "number");
      bonds_compiled += bond;
    } else if (type == "CUP") {
      check_bond(bond, n_sites, 1, false, "number");
      bonds_compiled += bond;
    } else if (type == "CDN") {
      check_bond(bond, n_sites, 1, false, "number");
      bonds_compiled += bond;
    } else {
      XDIAG_THROW(fmt::format("Invalid or undefined Bond type: \"{}\"", type));
    }
  }
  return bonds_compiled;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return BondList();
}

} // namespace xdiag::tj
