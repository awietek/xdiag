#include "compile.hpp"

#include <xdiag/operators/compiler.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/print_macro.hpp>

namespace xdiag::tj {

BondList compile(BondList const &bonds, double precision) try {
  using namespace operators;
  BondList bonds_explicit = compile_explicit(bonds, precision, "keep");
  BondList bonds_compiled;
  for (auto bond : bonds_explicit) {
    if (bond.type_defined()) {
      std::string type = bond.type();
      if (std::find(tj::special_bond_types.begin(),
                    tj::special_bond_types.end(),
                    type) == tj::special_bond_types.end()) {
        XDIAG_THROW(std::string("Invalid or undefined type: \"") + type +
                    std::string("\""));
      } else {

        // Exchange and Ising terms
        if ((type == "HB") || (type == "HEISENBERG")) {
          check_bond_has_correct_number_of_sites(bond, 2);
          check_bond_has_disjoint_sites(bond);
          bonds_compiled << Bond("ISING", bond.coupling(), bond.sites());
          bonds_compiled << Bond("EXCHANGE", bond.coupling(), bond.sites());
        } else if ((type == "TJHB") || (type == "TJHEISENBERG")) {
          check_bond_has_correct_number_of_sites(bond, 2);
          check_bond_has_disjoint_sites(bond);
          bonds_compiled << Bond("TJISING", bond.coupling(), bond.sites());
          bonds_compiled << Bond("EXCHANGE", bond.coupling(), bond.sites());
        } else if (type == "ISING") {
          check_bond_has_correct_number_of_sites(bond, 2);
          check_bond_has_disjoint_sites(bond);
          bonds_compiled << Bond("ISING", bond.coupling(), bond.sites());
        } else if (type == "TJISING") {
          check_bond_has_correct_number_of_sites(bond, 2);
          check_bond_has_disjoint_sites(bond);
          bonds_compiled << Bond("TJISING", bond.coupling(), bond.sites());
        } else if (type == "EXCHANGE") {
          check_bond_has_correct_number_of_sites(bond, 2);
          check_bond_has_disjoint_sites(bond);
          bonds_compiled << Bond("EXCHANGE", bond.coupling(), bond.sites());
        }

        // Hopping terms
        else if (type == "HOP") {
          check_bond_has_correct_number_of_sites(bond, 2);
          check_bond_has_disjoint_sites(bond);
          bonds_compiled << Bond("HOPUP", bond.coupling(), bond.sites());
          bonds_compiled << Bond("HOPDN", bond.coupling(), bond.sites());
        } else if (type == "HOPUP") {
          check_bond_has_correct_number_of_sites(bond, 2);
          check_bond_has_disjoint_sites(bond);
          bonds_compiled << Bond("HOPUP", bond.coupling(), bond.sites());
        } else if (type == "HOPDN") {
          check_bond_has_correct_number_of_sites(bond, 2);
          check_bond_has_disjoint_sites(bond);
          bonds_compiled << Bond("HOPDN", bond.coupling(), bond.sites());
        }

        // Number operators
        else if (type == "NUMBER") {
          check_bond_has_correct_number_of_sites(bond, 1);
          bonds_compiled << Bond("NUMBERUP", bond.coupling(), bond.sites());
          bonds_compiled << Bond("NUMBERDN", bond.coupling(), bond.sites());
        } else if (type == "NUMBERUP") {
          check_bond_has_correct_number_of_sites(bond, 1);
          bonds_compiled << Bond("NUMBERUP", bond.coupling(), bond.sites());
        } else if (type == "NUMBERDN") {
          check_bond_has_correct_number_of_sites(bond, 1);
          bonds_compiled << Bond("NUMBERDN", bond.coupling(), bond.sites());
        } else if (type == "SZ") {
          check_bond_has_correct_number_of_sites(bond, 1);
          bonds_compiled << Bond("NUMBERUP", 0.5 * bond.coupling(),
                                 bond.sites());
          bonds_compiled << Bond("NUMBERDN", -0.5 * bond.coupling(),
                                 bond.sites());
        }

        // Raising Lowering operators
        else if (type == "CDAGUP") {
          check_bond_has_correct_number_of_sites(bond, 1);
          bonds_compiled << Bond("CDAGUP", bond.coupling(), bond.sites());
        } else if (type == "CDAGDN") {
          check_bond_has_correct_number_of_sites(bond, 1);
          bonds_compiled << Bond("CDAGDN", bond.coupling(), bond.sites());
        } else if (type == "CUP") {
          check_bond_has_correct_number_of_sites(bond, 1);
          bonds_compiled << Bond("CUP", bond.coupling(), bond.sites());
        } else if (type == "CDN") {
          check_bond_has_correct_number_of_sites(bond, 1);
          bonds_compiled << Bond("CDN", bond.coupling(), bond.sites());
        }
      }
    }
  }
  return bonds_compiled;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return BondList();
}

} // namespace xdiag::tj
