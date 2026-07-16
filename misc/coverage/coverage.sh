#!/bin/bash
# SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
#
# SPDX-License-Identifier: Apache-2.0
#
# Build the xdiag tests with source-based coverage instrumentation, run them,
# and produce a terminal report plus a browsable HTML report.
#
# Usage:
#   misc/coverage/coverage.sh [catch2-test-filter]
#
# Examples:
#   misc/coverage/coverage.sh                 # full test suite
#   misc/coverage/coverage.sh "[linalg]"      # only tests tagged [linalg]
#
# See misc/coverage/README.md for details.

set -euo pipefail

# --- locate paths (repo root is two levels up from this script) --------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
BUILD_DIR="${REPO_ROOT}/build-cov"
HTML_DIR="${SCRIPT_DIR}/coverage-html"

TEST_FILTER="${1:-}"
NPROC=`nproc --all` 

# --- pick the toolchain's llvm-cov / llvm-profdata ---------------------------
# On macOS these must match Apple clang, so prefer the Xcode toolchain via
# xcrun; fall back to plain names elsewhere (e.g. Linux with clang).
if command -v xcrun >/dev/null 2>&1; then
  LLVM_PROFDATA=(xcrun llvm-profdata)
  LLVM_COV=(xcrun llvm-cov)
else
  LLVM_PROFDATA=(llvm-profdata)
  LLVM_COV=(llvm-cov)
fi

# Ignore third-party code and the test sources themselves in the report.
IGNORE_REGEX='(extern|tests|build|build-cov)/'

echo "==> Configuring coverage build in ${BUILD_DIR}"
cmake -B "${BUILD_DIR}" -S "${REPO_ROOT}" \
  -DCMAKE_BUILD_TYPE=Debug \
  -DXDIAG_COVERAGE=ON \
  -DBUILD_TESTING=ON

echo "==> Building tests"
cmake --build "${BUILD_DIR}" --target tests -j$NPROC

TEST_BIN="${BUILD_DIR}/tests/tests"
PROFRAW="${BUILD_DIR}/tests.profraw"
PROFDATA="${BUILD_DIR}/tests.profdata"

echo "==> Running tests${TEST_FILTER:+ (filter: ${TEST_FILTER})}"
LLVM_PROFILE_FILE="${PROFRAW}" "${TEST_BIN}" ${TEST_FILTER:+"${TEST_FILTER}"}

echo "==> Merging profile data"
"${LLVM_PROFDATA[@]}" merge -sparse "${PROFRAW}" -o "${PROFDATA}"

echo "==> Coverage report least-covered first"
"${LLVM_COV[@]}" report "${TEST_BIN}" \
  -instr-profile="${PROFDATA}" \
  -ignore-filename-regex="${IGNORE_REGEX}"

echo "==> Writing HTML report to ${HTML_DIR}"
"${LLVM_COV[@]}" show "${TEST_BIN}" \
  -instr-profile="${PROFDATA}" \
  -format=html -output-dir="${HTML_DIR}" \
  -ignore-filename-regex="${IGNORE_REGEX}" \
  -show-instantiation-summary

echo
echo "Done. Open the HTML report with:"
echo "  open ${HTML_DIR}/index.html"
