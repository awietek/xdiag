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
      ops_compiled += cpl * Op("SZSZ", op.sites());
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
      ops_compiled += cpl * Op("SZSZ", op.sites());
      ops_compiled += cpl * Op("EXCHANGE", op.sites());
    } else if (type == "TJSDOTS") {
      ops_compiled += cpl * Op("TJSZSZ", op.sites());
      ops_compiled += cpl * Op("EXCHANGE", op.sites());
    } else if (type == "HOP") {
      ops_compiled += cpl * Op("HOPUP", op.sites());
      ops_compiled += cpl * Op("HOPDN", op.sites());
    } else if (type == "NTOT") {
      ops_compiled += cpl * Op("NUP", op.sites());
      ops_compiled += cpl * Op("NDN", op.sites());
    } else if (type == "SZ") {
      ops_compiled += (Scalar(0.5) * cpl.scalar()) * Op("NUP", op.sites());
      ops_compiled +=
          (Scalar(-0.5) * cpl.scalar()) * Op("NDN", op.sites());
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
      ops_compiled += cpl * Op("SZSZ", op.sites());
      ops_compiled += cpl * Op("EXCHANGE", op.sites());
    } else if (type == "TJSDOTS") {
      ops_compiled += cpl * Op("TJSZSZ", op.sites());
      ops_compiled += cpl * Op("EXCHANGE", op.sites());
    } else if (type == "HOP") {
      ops_compiled += cpl * Op("HOPUP", op.sites());
      ops_compiled += cpl * Op("HOPDN", op.sites());
    } else if (type == "NTOT") {
      ops_compiled += cpl * Op("NUP", op.sites());
      ops_compiled += cpl * Op("NDN", op.sites());
    } else if (type == "SZ") {
      ops_compiled += (Scalar(0.5) * cpl.scalar()) * Op("NUP", op.sites());
      ops_compiled +=
          (Scalar(-0.5) * cpl.scalar()) * Op("NDN", op.sites());
    } else {
      ops_compiled += cpl * op;
    }
  }
  return ops_compiled;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

} // namespace xdiag::operators
