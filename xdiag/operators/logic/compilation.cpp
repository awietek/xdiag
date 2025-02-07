#include "compilation.hpp"

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
    if (type == "SdotS") {
      ops_compiled += cpl * Op("SzSz", op.sites());
      ops_compiled += cpl * Op("Exchange", op.sites());
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

    if (type == "SdotS") {
      ops_compiled += cpl * Op("SzSz", op.sites());
      ops_compiled += cpl * Op("Exchange", op.sites());
    } else if (type == "tJSdotS") {
      ops_compiled += cpl * Op("tJSzSz", op.sites());
      ops_compiled += cpl * Op("Exchange", op.sites());
    } else if (type == "Hop") {
      ops_compiled += cpl * Op("Hopup", op.sites());
      ops_compiled += cpl * Op("Hopdn", op.sites());
    } else if (type == "Ntot") {
      ops_compiled += cpl * Op("Nup", op.sites());
      ops_compiled += cpl * Op("Ndn", op.sites());
    } else if (type == "Sz") {
      ops_compiled += (Scalar(0.5) * cpl.scalar()) * Op("Nup", op.sites());
      ops_compiled += (Scalar(-0.5) * cpl.scalar()) * Op("Ndn", op.sites());
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

    if (type == "SdotS") {
      ops_compiled += cpl * Op("SzSz", op.sites());
      ops_compiled += cpl * Op("Exchange", op.sites());
    } else if (type == "tJSdotS") {
      ops_compiled += cpl * Op("tJSzSz", op.sites());
      ops_compiled += cpl * Op("Exchange", op.sites());
    } else if (type == "Hop") {
      ops_compiled += cpl * Op("Hopup", op.sites());
      ops_compiled += cpl * Op("Hopdn", op.sites());
    } else if (type == "Ntot") {
      ops_compiled += cpl * Op("Nup", op.sites());
      ops_compiled += cpl * Op("Ndn", op.sites());
    } else if (type == "Sz") {
      ops_compiled += (Scalar(0.5) * cpl.scalar()) * Op("Nup", op.sites());
      ops_compiled += (Scalar(-0.5) * cpl.scalar()) * Op("Ndn", op.sites());
    } else {
      ops_compiled += cpl * op;
    }
  }
  return ops_compiled;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

template <> OpSum compile<Spinhalf>(OpSum const &ops) try {
  return compile_spinhalf(ops);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}
template <> OpSum compile<tJ>(OpSum const &ops) try {
  return compile_tj(ops);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}
template <> OpSum compile<Electron>(OpSum const &ops) try {
  return compile_electron(ops);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}
#ifdef XDIAG_USE_MPI
template <> OpSum compile<SpinhalfDistributed>(OpSum const &ops) try {
  return compile_spinhalf(ops);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}
template <> OpSum compile<tJDistributed>(OpSum const &ops) try {
  return compile_tj(ops);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}
#endif

} // namespace xdiag::operators
