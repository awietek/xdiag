#include "compiler.hpp"

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

} // namespace xdiag::operators
