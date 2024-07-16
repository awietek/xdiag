#include "compile.hpp"

#include <xdiag/operators/compiler.hpp>
#include <xdiag/operators/non_branching_op.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/print_macro.hpp>
#include <xdiag/utils/timing.hpp>

namespace xdiag::spinhalf {

OpSum compile(OpSum const &ops, int64_t n_sites, double precision) try {
  using namespace operators;
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

} // namespace xdiag::spinhalf
