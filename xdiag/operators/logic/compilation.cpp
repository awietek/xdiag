#include "compilation.hpp"
#include <xdiag/operators/logic/valid.hpp>
#include <xdiag/utils/scalar.hpp>

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

  // Additional changes for Spinhalf
  OpSum ops_double;
  for (auto const &[cpl, op] : ops_compiled) {
    std::string type = op.type();

    // Replacement if of is on same site
    if ((type == "SzSz") && (op[0] == op[1])) {
      auto cpl2 = Scalar(0.25) * cpl.scalar();
      ops_double += cpl2 * Op("Id");
    } else if ((type == "Exchange") && (op[0] == op[1])) {
      auto cpl2 = Scalar(0.5) * cpl.scalar();
      ops_double += cpl2 * Op("Id");
    } else if ((type == "SdotS") && (op[0] == op[1])) {
      auto cpl2 = Scalar(0.75) * cpl.scalar();
      ops_double += cpl2 * Op("Id");
    } else {
      ops_double += cpl * op;
    }
  }

  // Compbine Matrix operators on same sites
  OpSum ops_matrix;
  std::map<std::vector<int64_t>, Matrix> matrix_on_sites;
  for (auto const &[cpl, op] : ops_double) {
    if (op.type() == "Matrix") {
      auto sites = op.sites();
      auto coeff = cpl.scalar();
      if (matrix_on_sites.count(sites)) { // sites already exist
        matrix_on_sites[sites] += op.matrix() * coeff;
      } else {
        matrix_on_sites[sites] = op.matrix() * coeff;
      }
    } else {
      ops_matrix += cpl * op;
    }
  }

  for (auto const &[sites, matrix] : matrix_on_sites) {
    ops_matrix += Op("Matrix", sites, matrix);
  }

  return ops_matrix;
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

  // Additional changes for tJ
  OpSum ops_final;
  for (auto const &[cpl, op] : ops_compiled) {
    std::string type = op.type();

    // Replacement if of is on same site
    if ((type == "SzSz") && (op[0] == op[1])) {
      auto cpl2 = Scalar(0.25) * cpl.scalar();
      ops_final += cpl2 * Op("Nup", op[0]);
      ops_final += cpl2 * Op("Ndn", op[0]);
    } else if ((type == "tJSzSz") && (op[0] == op[1])) {
      // do nothing, this operator is identical 0
    } else if ((type == "Exchange") && (op[0] == op[1])) {
      auto cpl2 = Scalar(0.5) * cpl.scalar();
      ops_final += cpl2 * Op("Nup", op[0]);
      ops_final += cpl2 * Op("Ndn", op[0]);
    } else if ((type == "SdotS") && (op[0] == op[1])) {
      auto cpl2 = Scalar(0.75) * cpl.scalar();
      ops_final += cpl2 * Op("Nup", op[0]);
      ops_final += cpl2 * Op("Ndn", op[0]);
    } else if ((type == "tJSdotS") && (op[0] == op[1])) {
      auto cpl2 = Scalar(0.5) * cpl.scalar();
      ops_final += cpl2 * Op("Nup", op[0]);
      ops_final += cpl2 * Op("Ndn", op[0]);
    } else if ((type == "NtotNtot") && (op[0] == op[1])) {
      ops_final += cpl * Op("Nup", op[0]);
      ops_final += cpl * Op("Ndn", op[0]);
    } else if ((type == "Hopup") && (op[0] == op[1])) {
      auto cpl2 = Scalar(-2.0) * cpl.scalar();
      ops_final += cpl2 * Op("Nup", op[0]);
    } else if ((type == "Hopdn") && (op[0] == op[1])) {
      auto cpl2 = Scalar(-2.0) * cpl.scalar();
      ops_final += cpl2 * Op("Ndn", op[0]);
    } else {
      ops_final += cpl * op;
    }
  }
  return ops_final;
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

  // Additional changes for Electron
  OpSum ops_final;
  for (auto const &[cpl, op] : ops_compiled) {
    std::string type = op.type();
    auto s = cpl.scalar();
    // Replacement if of is on same site
    if ((type == "SzSz") && (op[0] == op[1])) {
      ops_final += (Scalar(0.25) * s) * Op("Nup", op[0]);
      ops_final += (Scalar(0.25) * s) * Op("Ndn", op[0]);
      ops_final += (Scalar(-0.5) * s) * Op("Nupdn", op[0]);
    } else if ((type == "Exchange") && (op[0] == op[1])) {
      ops_final += (Scalar(0.5) * s) * Op("Nup", op[0]);
      ops_final += (Scalar(0.5) * s) * Op("Ndn", op[0]);
      ops_final += (Scalar(-1.0) * s) * Op("Nupdn", op[0]);
    } else if ((type == "SdotS") && (op[0] == op[1])) {
      ops_final += (Scalar(0.5) * s) * Op("Nup", op[0]);
      ops_final += (Scalar(0.5) * s) * Op("Ndn", op[0]);
      ops_final += (Scalar(-1.0) * s) * Op("Nupdn", op[0]);
    } else if ((type == "NtotNtot") && (op[0] == op[1])) {
      ops_final += cpl * Op("Nup", op[0]);
      ops_final += cpl * Op("Ndn", op[0]);
      ops_final += (Scalar(2.0) * cpl.scalar()) * Op("Nupdn", op[0]);
    } else if ((type == "NupdnNupdn") && (op[0] == op[1])) {
      ops_final += cpl * Op("Nupdn", op[0]);
    } else if ((type == "Hopup") && (op[0] == op[1])) {
      auto cpl2 = Scalar(-2.0) * cpl.scalar();
      ops_final += cpl2 * Op("Nup", op[0]);
    } else if ((type == "Hopdn") && (op[0] == op[1])) {
      auto cpl2 = Scalar(-2.0) * cpl.scalar();
      ops_final += cpl2 * Op("Ndn", op[0]);
    } else {
      ops_final += cpl * op;
    }
  }
  return ops_final;
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
