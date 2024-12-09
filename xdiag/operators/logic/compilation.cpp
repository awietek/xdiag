#include "compilation.hpp"

#include <xdiag/operators/non_branching_op.hpp>

namespace xdiag::operators {

OpSum clean_zeros(OpSum const &ops) try {
  OpSum ops_clean;
  for (auto const &[cpl, op] : ops.plain()) {
    if (cpl.scalar() != zero(cpl.scalar())) {
      ops_clean += cpl * op;
    }
  }
  return ops_clean;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

OpSum compile_spinhalf(OpSum const &ops) try {
  OpSum ops_clean = clean_zeros(ops.plain());
  OpSum ops_compiled;
  for (auto [cpl, op] : ops_clean) {
    std::string type = op.type();
    if (type == "SDOTS") {
      ops_compiled += cpl * Op("ISING", op.sites());
      ops_compiled += cpl * Op("EXCHANGE", op.sites());
    } else {
      ops_compiled += cpl * op;
    }
  }
  return ops_compiled;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

OpSum compile_tj(OpSum const &ops) try {
  OpSum ops_clean = clean_zeros(ops);
  OpSum ops_compiled;

  for (auto [cpl, op] : ops_clean) {
    std::string type = op.type();

    if (type == "SDOTS") {
      ops_compiled += cpl * Op("ISING", op.sites());
      ops_compiled += cpl * Op("EXCHANGE", op.sites());
    } else if (type == "TJSDOTS") {
      ops_compiled += cpl * Op("TJISING", op.sites());
      ops_compiled += cpl * Op("EXCHANGE", op.sites());
    } else if (type == "HOP") {
      ops_compiled += cpl * Op("HOPUP", op.sites());
      ops_compiled += cpl * Op("HOPDN", op.sites());
    } else if (type == "NUMBER") {
      ops_compiled += cpl * Op("NUMBERUP", op.sites());
      ops_compiled += cpl * Op("NUMBERDN", op.sites());
    } else if (type == "SZ") {
      ops_compiled += (Scalar(0.5) * cpl.scalar()) * Op("NUMBERUP", op.sites());
      ops_compiled +=
          (Scalar(-0.5) * cpl.scalar()) * Op("NUMBERDN", op.sites());
    } else {
      ops_compiled += cpl * op;
    }
  }
  return ops_compiled;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

OpSum compile_electron(OpSum const &ops) try {
  OpSum ops_clean = clean_zeros(ops);
  OpSum ops_compiled;

  for (auto [cpl, op] : ops_clean) {
    std::string type = op.type();

    if (type == "SDOTS") {
      ops_compiled += cpl * Op("ISING", op.sites());
      ops_compiled += cpl * Op("EXCHANGE", op.sites());
    } else if (type == "TJSDOTS") {
      ops_compiled += cpl * Op("TJISING", op.sites());
      ops_compiled += cpl * Op("EXCHANGE", op.sites());
    } else if (type == "HOP") {
      ops_compiled += cpl * Op("HOPUP", op.sites());
      ops_compiled += cpl * Op("HOPDN", op.sites());
    } else if (type == "NUMBER") {
      ops_compiled += cpl * Op("NUMBERUP", op.sites());
      ops_compiled += cpl * Op("NUMBERDN", op.sites());
    } else if (type == "SZ") {
      ops_compiled += (Scalar(0.5) * cpl.scalar()) * Op("NUMBERUP", op.sites());
      ops_compiled +=
          (Scalar(-0.5) * cpl.scalar()) * Op("NUMBERDN", op.sites());
    } else {
      ops_compiled += cpl * op;
    }
  }
  return ops_compiled;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

} // namespace xdiag::operators
