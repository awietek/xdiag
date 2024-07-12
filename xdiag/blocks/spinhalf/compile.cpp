#include "compile.hpp"

#include <xdiag/operators/compiler.hpp>
#include <xdiag/operators/non_branching_bonds.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/print_macro.hpp>
#include <xdiag/utils/timing.hpp>

namespace xdiag::spinhalf {

BondList compile(BondList const &bonds, int64_t n_sites, double precision) try {
  using namespace operators;

  BondList bonds_explicit = make_explicit(bonds);
  BondList bonds_compiled;
  for (auto bond : bonds_explicit) {

    if (bond.ismatrix()) {
      BondList bonds_nb = non_branching_bonds(bond, precision);
      bonds_compiled += bonds_nb;
    } else {
      std::string type = bond.type();

      if (type == "HB") {
        check_bond(bond, n_sites, 2, true, "number");
        bonds_compiled += Bond("ISING", bond.coupling(), bond.sites());
        bonds_compiled += Bond("EXCHANGE", bond.coupling(), bond.sites());
      } else if (type == "ISING") {
        check_bond(bond, n_sites, 2, true, "number");
        bonds_compiled += bond;
      } else if (type == "EXCHANGE") {
        check_bond(bond, n_sites, 2, true, "number");
        bonds_compiled += bond;
      } else if (type == "SZ") {
        check_bond(bond, n_sites, 1, false, "number");
        bonds_compiled += bond;
      } else if (type == "S+") {
        check_bond(bond, n_sites, 1, false, "number");
        bonds_compiled += bond;
      } else if (type == "S-") {
        check_bond(bond, n_sites, 1, false, "number");
        bonds_compiled += bond;
      } else if (type == "SCALARCHIRALITY") {
        check_bond(bond, n_sites, 3, true, "number");
        bonds_compiled += bond;
      } else {
        XDIAG_THROW(fmt::format("Invalid or undefined type: \"{}\"", type));
      }
    }
  }
  return bonds_compiled;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return BondList();
}

} // namespace xdiag::spinhalf
