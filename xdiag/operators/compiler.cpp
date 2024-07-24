#include "compiler.hpp"

#include <xdiag/operators/non_branching_op.hpp>

namespace xdiag::operators {

OpSum clean_zeros(OpSum const &ops, double precision) try {
  OpSum clean_ops;
  for (auto const &op : ops) {

    // if Op is only defined by string, cannot be cleaned
    if (!op.isexplicit()) {
      XDIAG_THROW("Cannot clean zero Op since is is not explicit. Maybe call "
                  "\"make_explicit\" first");
    }

    Coupling cpl = op.coupling();
    if (cpl.is<double>()) {
      double cval = cpl.as<double>();
      if (std::abs(cval) > precision) {
        clean_ops += op;
      }
    } else if (cpl.is<complex>()) {
      complex cval = cpl.as<complex>();
      if (std::abs(cval) > precision) {
        clean_ops += op;
      }
    }
    // coupling is matrix
    else {
      clean_ops += op;
    }
  }
  return clean_ops;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void check_op(Op const &op, int64_t n_sites_total, int64_t n_sites_op,
              bool disjoint, std::string type) {
  check_op_in_range(op, n_sites_total);
  check_op_has_correct_number_of_sites(op, n_sites_op);
  if (disjoint) {
    check_op_has_disjoint_sites(op);
  }
  if (type == "number") {
    check_op_coupling_has_type(op, "double", "complex");
  } else if (type == "matrix") {
    check_op_coupling_has_type(op, "arma::mat", "arma::cx_mat");
  } else {
    check_op_coupling_has_type(op, type);
  }
}

void check_op_in_range(Op const &op, int64_t n_sites) try {
  for (int64_t s : op.sites()) {
    if (s >= n_sites) {
      XDIAG_THROW(
          fmt::format("Op site number {} exceeds range of number of sites = {}",
                      s, n_sites));
    } else if (s < 0) {
      XDIAG_THROW(fmt::format("Op site number {} found to be < 0", s));
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void check_ops_in_range(OpSum const &ops, int64_t n_sites) try {
  for (Op op : ops) {
    check_op_in_range(op, n_sites);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void check_op_has_correct_number_of_sites(Op const &op, int64_t ns) try {
  if (op.size() != ns) {
    std::string msg = std::string("Op of type \"") + op.type() +
                      std::string("\": number of sites given is ") +
                      std::to_string(op.size()) +
                      std::string(" while expected ") + std::to_string(ns);
    XDIAG_THROW(msg);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void check_op_has_disjoint_sites(Op const &op) try {
  if (!sites_disjoint(op)) {
    XDIAG_THROW(std::string("Op of type \"") + op.type() +
                std::string("\": sites are not disjoint as expected."));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void check_op_coupling_has_type(Op const &op, std::string type) try {
  std::string op_type = op.coupling().type();
  if (!(op_type == type)) {
    XDIAG_THROW(fmt::format("Coupling of Op is expected to be of type \"{}\" "
                            "but received a type \"{}\"",
                            type, op_type));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void check_op_coupling_has_type(Op const &op, std::string type1,
                                std::string type2) try {
  std::string op_type = op.coupling().type();
  if (!(op_type == type1) && !(op_type == type2)) {
    XDIAG_THROW(fmt::format("Coupling of Op is expected to be of type either "
                            "\"{}\" or \"{}\" but received a type \"{}\"",
                            type1, type2, op_type));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

OpSum compile(OpSum const &ops, Spinhalf const &block, double precision) try {
  int64_t n_sites = block.n_sites();
  OpSum ops_explicit = make_explicit(ops);
  OpSum ops_clean = clean_zeros(ops_explicit, precision);

  OpSum ops_compiled;
  for (auto op : ops_clean) {

    if (op.ismatrix()) {
      OpSum ops_nb = non_branching_ops(op, precision);
      ops_compiled += ops_nb;
    } else {
      std::string type = op.type();

      if (type == "HB") {
        check_op(op, n_sites, 2, true, "number");
        ops_compiled += Op("ISING", op.coupling(), op.sites());
        ops_compiled += Op("EXCHANGE", op.coupling(), op.sites());
      } else if (type == "ISING") {
        check_op(op, n_sites, 2, true, "number");
        ops_compiled += op;
      } else if (type == "EXCHANGE") {
        check_op(op, n_sites, 2, true, "number");
        ops_compiled += op;
      } else if (type == "SZ") {
        check_op(op, n_sites, 1, false, "number");
        ops_compiled += op;
      } else if (type == "S+") {
        check_op(op, n_sites, 1, false, "number");
        ops_compiled += op;
      } else if (type == "S-") {
        check_op(op, n_sites, 1, false, "number");
        ops_compiled += op;
      } else if (type == "SCALARCHIRALITY") {
        check_op(op, n_sites, 3, true, "number");
        ops_compiled += op;
      } else {
        XDIAG_THROW(fmt::format("Invalid or undefined type: \"{}\"", type));
      }
    }
  }
  return ops_compiled;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

OpSum compile(OpSum const &ops, tJ const &block, double precision) try {
  int64_t n_sites = block.n_sites();
  OpSum ops_explicit = make_explicit(ops);
  OpSum ops_clean = clean_zeros(ops_explicit, precision);

  OpSum ops_compiled;
  for (auto op : ops_clean) {
    std::string type = op.type();

    // Exchange and Ising terms
    if (type == "HB") {
      check_op(op, n_sites, 2, true, "number");
      ops_compiled += Op("ISING", op.coupling(), op.sites());
      ops_compiled += Op("EXCHANGE", op.coupling(), op.sites());
    } else if (type == "TJHB") {
      check_op(op, n_sites, 2, true, "number");
      ops_compiled += Op("TJISING", op.coupling(), op.sites());
      ops_compiled += Op("EXCHANGE", op.coupling(), op.sites());
    } else if (type == "ISING") {
      check_op(op, n_sites, 2, true, "number");
      ops_compiled += op;
    } else if (type == "TJISING") {
      check_op(op, n_sites, 2, true, "number");
      ops_compiled += op;
    } else if (type == "EXCHANGE") {
      check_op(op, n_sites, 2, true, "number");
      ops_compiled += op;
    }

    // Hopping terms
    else if (type == "HOP") {
      check_op(op, n_sites, 2, true, "number");
      ops_compiled += Op("HOPUP", op.coupling(), op.sites());
      ops_compiled += Op("HOPDN", op.coupling(), op.sites());
    } else if (type == "HOPUP") {
      check_op(op, n_sites, 2, true, "number");
      ops_compiled += op;
    } else if (type == "HOPDN") {
      check_op(op, n_sites, 2, true, "number");
      ops_compiled += op;
    }

    // Number operators
    else if (type == "NUMBER") {
      check_op(op, n_sites, 1, false, "number");
      ops_compiled += Op("NUMBERUP", op.coupling(), op.sites());
      ops_compiled += Op("NUMBERDN", op.coupling(), op.sites());
    } else if (type == "NUMBERUP") {
      check_op(op, n_sites, 1, false, "number");
      ops_compiled += Op("NUMBERUP", op.coupling(), op.sites());
    } else if (type == "NUMBERDN") {
      check_op(op, n_sites, 1, false, "number");
      ops_compiled += Op("NUMBERDN", op.coupling(), op.sites());
    } else if (type == "SZ") {
      check_op(op, n_sites, 1, false, "number");
      if (op.coupling().is<double>()) {
        ops_compiled +=
            Op("NUMBERUP", 0.5 * op.coupling().as<double>(), op.sites());
        ops_compiled +=
            Op("NUMBERDN", -0.5 * op.coupling().as<double>(), op.sites());
      } else if (op.coupling().is<complex>()) {
        ops_compiled +=
            Op("NUMBERUP", 0.5 * op.coupling().as<complex>(), op.sites());
        ops_compiled +=
            Op("NUMBERDN", -0.5 * op.coupling().as<complex>(), op.sites());
      }
    }

    // Raising Lowering operators
    else if (type == "CDAGUP") {
      check_op(op, n_sites, 1, false, "number");
      ops_compiled += op;
    } else if (type == "CDAGDN") {
      check_op(op, n_sites, 1, false, "number");
      ops_compiled += op;
    } else if (type == "CUP") {
      check_op(op, n_sites, 1, false, "number");
      ops_compiled += op;
    } else if (type == "CDN") {
      check_op(op, n_sites, 1, false, "number");
      ops_compiled += op;
    } else {
      XDIAG_THROW(fmt::format("Invalid or undefined Op type: \"{}\"", type));
    }
  }
  return ops_compiled;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

OpSum compile(OpSum const &ops, Electron const &block, double precision) try {
  int64_t n_sites = block.n_sites();
  OpSum ops_explicit = make_explicit(ops);
  OpSum ops_clean = clean_zeros(ops_explicit, precision);

  OpSum ops_compiled;
  for (auto op : ops_clean) {
    std::string type = op.type();

    // Exchange and Ising terms
    if (type == "HB") {
      check_op(op, n_sites, 2, true, "number");
      ops_compiled += Op("ISING", op.coupling(), op.sites());
      ops_compiled += Op("EXCHANGE", op.coupling(), op.sites());
    } else if (type == "ISING") {
      check_op(op, n_sites, 2, true, "number");
      ops_compiled += Op("ISING", op.coupling(), op.sites());
    } else if (type == "EXCHANGE") {
      check_op(op, n_sites, 2, true, "number");
      ops_compiled += Op("EXCHANGE", op.coupling(), op.sites());

      // Hopping terms
    } else if (type == "HOP") {
      check_op(op, n_sites, 2, true, "number");
      ops_compiled += Op("HOPUP", op.coupling(), op.sites());
      ops_compiled += Op("HOPDN", op.coupling(), op.sites());
    } else if (type == "HOPUP") {
      check_op(op, n_sites, 2, true, "number");
      ops_compiled += Op("HOPUP", op.coupling(), op.sites());
    } else if (type == "HOPDN") {
      check_op(op, n_sites, 2, true, "number");
      ops_compiled += Op("HOPDN", op.coupling(), op.sites());
    }

    // Number operators
    else if (type == "NUMBER") {
      check_op(op, n_sites, 1, false, "number");
      ops_compiled += Op("NUMBERUP", op.coupling(), op.sites());
      ops_compiled += Op("NUMBERDN", op.coupling(), op.sites());
    } else if (type == "NUMBERUP") {
      check_op(op, n_sites, 1, false, "number");
      ops_compiled += Op("NUMBERUP", op.coupling(), op.sites());
    } else if (type == "NUMBERDN") {
      check_op(op, n_sites, 1, false, "number");
      ops_compiled += Op("NUMBERDN", op.coupling(), op.sites());
    } else if (type == "SZ") {
      check_op(op, n_sites, 1, false, "number");
      if (op.coupling().is<double>()) {
        ops_compiled +=
            Op("NUMBERUP", 0.5 * op.coupling().as<double>(), op.sites());
        ops_compiled +=
            Op("NUMBERDN", -0.5 * op.coupling().as<double>(), op.sites());
      } else if (op.coupling().is<complex>()) {
        ops_compiled +=
            Op("NUMBERUP", 0.5 * op.coupling().as<complex>(), op.sites());
        ops_compiled +=
            Op("NUMBERDN", -0.5 * op.coupling().as<complex>(), op.sites());
      }
    }

    // Raising Lowering operators
    else if (type == "CDAGUP") {
      check_op(op, n_sites, 1, false, "number");
      ops_compiled += Op("CDAGUP", op.coupling(), op.sites());
    } else if (type == "CDAGDN") {
      check_op(op, n_sites, 1, false, "number");
      ops_compiled += Op("CDAGDN", op.coupling(), op.sites());
    } else if (type == "CUP") {
      check_op(op, n_sites, 1, false, "number");
      ops_compiled += Op("CUP", op.coupling(), op.sites());
    } else if (type == "CDN") {
      check_op(op, n_sites, 1, false, "number");
      ops_compiled += Op("CDN", op.coupling(), op.sites());
    } else {
      XDIAG_THROW(fmt::format("Invalid or undefined Op type: \"{}\"", type));
    }
  }
  // Set Hubbbbard U term again
  if (ops.defined("U")) {
    ops_compiled["U"] = ops["U"];
  }

  return ops_compiled;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return OpSum();
}

#ifdef XDIAG_USE_MPI
OpSum compile(OpSum const &ops, SpinhalfDistributed const &block,
              double precision) {
  int64_t n_sites = block.n_sites();
  int64_t n_up = block.n_up();
  return compile(ops, Spinhalf(n_sites, n_up), precision);
}
OpSum compile(OpSum const &ops, tJDistributed const &block, double precision) {
  int64_t n_sites = block.n_sites();
  int64_t n_up = block.n_up();
  int64_t n_dn = block.n_dn();
  return compile(ops, tJ(n_sites, n_up, n_dn), precision);
}
#endif

} // namespace xdiag::operators
