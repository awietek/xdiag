// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "fermi_string.hpp"

#include <algorithm>
#include <xdiag/operators/logic/valid.hpp>
#include <xdiag/utils/split.hpp>

namespace xdiag::operators {

static constexpr std::vector<std::string> valid_tokens = {"dagup", "dagdn",
                                                          "up", "dn"};

int64_t check_fermi_string(Op const &op) try {
  check_valid(op);

  std::vector<std::string> tokens = utils::split(op.type(), "C");
  for (auto t : tokens) {
    if (std::find(valid_tokens.begin(), valid_tokens.end(), t) ==
        valid_tokens.end()) {
      XDIAG_THROW(fmt::format("Invalid token found: \"{}\"", t));
    }
  }
  if (!op.hassites()) {
    XDIAG_THROW("Op defining a FermiString needs sites defined.");
  } else if (op.hasmatrix()) {
    XDIAG_THROW("Op defining a FermiString cannot have a matrix defined.");
  } else if (op.sites().size() != tokens.size()) {
    XDIAG_THROW("Number of sites does not match the number of Fermi operators");
  }
  return tokens.size();
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool is_fermi_string(Op const &op) {
  std::vector<std::string> tokens = utils::split(op.type(), "C");
  for (auto t : tokens) {
    if (std::find(valid_tokens.begin(), valid_tokens.end(), t) ==
        valid_tokens.end()) {
      return false;
    }
  }
  return true;
}

template <typename bit_t>
FermiString<bit_t>::FermiString(Op const &op) try
    : size_(check_fermi_string(op)) operator_(size_), set_mask_(size_), {
  check_fermi_string(op);

  std::vector<std::string> tokens = utils::split(op.type(), "C");
  std::vector<int64_t> sites = op.sites();

  for (int64_t i = 0; i < size_; ++i) {
    std::string t = tokens[i];
    if (t == "dagup") {
      operator_[i] = cdagup;
    } else if (t == "dagdn") {
      operator_[i] = cdagdn;
    } else if (t == "up") {
      operator_[i] = cup;
    } else if (t == "dn") {
      operator_[i] = cdn;
    } else {
      XDIAG_THROW("Invalid token");
    }
    set_mask_[i] = (bit_t)1 << sites[i];
    fermi_mask_[i] = set_mask_[i] - 1;
  }

} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename bit_t>
bool FermiString<bit_t>::non_zero_term(bit_t up, bit_t dn) const {
  for (int64_t i = 0; i < size_; ++i) {
    switch (operator_[i]) {
    case cdagup:
      if (!(up & set_mask_[i])) {
        return false;
      }
      break;
    case cdagdn:
      if (!(dn & set_mask_[i])) {
        return false;
      }
      break;
    case cup:
      if (up & set_mask_[i]) {
        return false;
      }
      break;
    case cdn:
      if (dn & set_mask_[i]) {
        return false;
      }
      break;
    }
  }
  return true;
}

} // namespace xdiag::operators
