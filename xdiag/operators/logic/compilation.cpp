// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "compilation.hpp"
#include <xdiag/operators/logic/valid.hpp>
#include <xdiag/operators/logic/order.hpp>
#include <xdiag/utils/scalar.hpp>

#include <xdiag/blocks/electron.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/blocks/tj.hpp>

#ifdef XDIAG_USE_MPI
#include <xdiag/blocks/electron_distributed.hpp>
#include <xdiag/blocks/spinhalf_distributed.hpp>
#include <xdiag/blocks/tj_distributed.hpp>
#endif

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
  OpSum ops_clean = clean_zeros(order(ops));
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

  return ops_double;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

OpSum compile_tj(OpSum const &ops) try {
  OpSum ops_clean = clean_zeros(order(ops));
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
  OpSum ops_clean = clean_zeros(order(ops));
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
    } else if (type == "Nhup") {
      ops_compiled += cpl * Op("Id");
      ops_compiled -= cpl * Op("Nup", op.sites());
    } else if (type == "Nhdn") {
      ops_compiled += Op("Id");
      ops_compiled -= cpl * Op("Ndn", op.sites());
    } else if (type == "Nhtot") {
      ops_compiled += (Scalar(2.0) * cpl.scalar()) * Op("Id");
      ops_compiled -= cpl * Op("Nup", op.sites());
      ops_compiled -= cpl * Op("Ndn", op.sites());
    } else if (type == "NhupNdn") {
      ops_compiled += cpl * Op("Ndn", op[1]);
      ops_compiled -= cpl * Op("NupNdn", op.sites());
    } else if (type == "NupNhdn") {
      ops_compiled += cpl * Op("Nup", op[0]);
      ops_compiled -= cpl * Op("NupNdn", op.sites());
    } else if (type == "NhupNhdn") {
      ops_compiled += cpl * Op("Id");
      ops_compiled -= cpl * Op("Nup", op[0]);
      ops_compiled -= cpl * Op("Ndn", op[1]);
      ops_compiled += cpl * Op("NupNdn", op.sites());
    } else if (type == "NhupNhup") {
      // Special case
      if (op[0] == op[1]) {
          ops_compiled += cpl * Op("Id");
          ops_compiled -= cpl * Op("Nup", op[0]);
          continue;
      }
      ops_compiled += cpl * Op("Id");
      ops_compiled -= cpl * Op("Nup", op[0]);
      ops_compiled -= cpl * Op("Nup", op[1]);
      ops_compiled += cpl * Op("NupNup", op.sites());
    } else if (type == "NhdnNhdn") {
      // Special case
      if (op[0] == op[1]) {
          ops_compiled += cpl * Op("Id");
          ops_compiled -= cpl * Op("Ndn", op[0]);
          continue;
      }
      ops_compiled += cpl * Op("Id");
      ops_compiled -= cpl * Op("Ndn", op[0]);
      ops_compiled -= cpl * Op("Ndn", op[1]);
      ops_compiled += cpl * Op("NdnNdn", op.sites());
    } else if (type == "NhupNup") {
      // Special case
      if (op[0] == op[1]) {
          continue;
      }
      ops_compiled += cpl * Op("Nup", op[1]);
      ops_compiled -= cpl * Op("NupNup", op.sites());
    } else if (type == "NhdnNdn") {
      // Special case
      if (op[0] == op[1]) {
          continue;
      }
      ops_compiled += cpl * Op("Ndn", op[1]);
      ops_compiled -= cpl * Op("NdnNdn", op.sites());
    } else if ((type == "NhtotNhtot")) {
      auto s2 = Scalar(2.0) * cpl.scalar();
      ops_compiled += (s2 + Scalar(2.0)) * Op("Id");
      ops_compiled -= s2 * Op("Nup", op[0]);
      ops_compiled -= s2 * Op("Ndn", op[0]);
      ops_compiled -= s2 * Op("Nup", op[1]);
      ops_compiled -= s2 * Op("Ndn", op[1]);
      ops_compiled += cpl * Op("NtotNtot", op.sites());
    } else if ((type == "NhtotNtot") && (op[0] == op[1])) {
      ops_compiled += cpl * Op("Nup", op[0]);
      ops_compiled += cpl * Op("Ndn", op[0]);
      ops_compiled -= (Scalar(2.0) * cpl.scalar()) * Op("Nupdn", op[0]);
    } else if (type == "NhtotNtot") {
      auto s2 = Scalar(2.0) * cpl.scalar();
      ops_compiled += s2 * Op("Nup", op[1]);
      ops_compiled += s2 * Op("Ndn", op[1]);
      ops_compiled -= cpl * Op("NtotNtot", op.sites());
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
    } else if ((type == "NupNdn") && (op[0] == op[1])) {
        ops_final += cpl * Op("Nupdn", op[0]);
    } else if ((type == "NupNup") && (op[0] == op[1])) {
        ops_final += cpl * Op("Nup", op[0]);
    } else if ((type == "NdnNdn") && (op[0] == op[1])) {
        ops_final += cpl * Op("Ndn", op[0]);
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
template <> OpSum compile<ElectronDistributed>(OpSum const &ops) try {
  return compile_electron(ops);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}
#endif

} // namespace xdiag::operators
