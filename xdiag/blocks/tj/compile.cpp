#include "compile.hpp"

#include <xdiag/operators/compiler.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/print_macro.hpp>

namespace xdiag::tj {

OpSum compile(OpSum const &ops, int64_t n_sites, double precision) try {
  using namespace operators;
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
        ops_compiled += Op("NUMBERDN", -0.5 * op.coupling().as<complex>(),
                               op.sites());
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
  return OpSum();
}

} // namespace xdiag::tj
