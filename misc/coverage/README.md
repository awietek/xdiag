<!--
SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>

SPDX-License-Identifier: Apache-2.0
-->

# Test coverage

xdiag uses **source-based coverage** (LLVM's `-fprofile-instr-generate
-fcoverage-mapping`), which is the right choice for the Clang toolchains used to
build the library (Apple Clang on macOS, upstream Clang elsewhere). This gives
accurate line-, region-, and branch-level coverage for templated and header-heavy
code.

## Quick start

From the repository root:

```bash
misc/coverage/coverage.sh
```

This will:

1. Configure a separate build directory `build-cov/` with `-DXDIAG_COVERAGE=ON`
   (your normal `build/` is left untouched).
2. Build and run the `tests` executable, emitting raw profile data.
3. Print a terminal report, least-covered files first.
4. Write a browsable HTML report to `misc/coverage/coverage-html/`
   (gitignored — see below).

Open the HTML report with:

```bash
open misc/coverage/coverage-html/index.html
```

### Covering just one area

Pass a [Catch2](https://github.com/catchorg/Catch2) test filter to restrict the
run (tag or test-name pattern):

```bash
misc/coverage/coverage.sh "[linalg]"
misc/coverage/coverage.sh "lobpcg*"
```

The report then reflects only the code exercised by the selected tests.

## How it works / doing it by hand

The `XDIAG_COVERAGE` CMake option (see the top-level `CMakeLists.txt`) attaches
the instrumentation flags to the `xdiag` library target with `PUBLIC` visibility,
so the `tests` executable — which links the library — inherits them
automatically. The equivalent manual steps are:

```bash
# 1. Configure + build with instrumentation
cmake -B build-cov -DCMAKE_BUILD_TYPE=Debug -DXDIAG_COVERAGE=ON -DBUILD_TESTING=ON
cmake --build build-cov --target tests -j

# 2. Run tests, emitting a raw profile
LLVM_PROFILE_FILE=build-cov/tests.profraw ./build-cov/tests/tests

# 3. Merge the raw profile
xcrun llvm-profdata merge -sparse build-cov/tests.profraw -o build-cov/tests.profdata

# 4a. Terminal summary
xcrun llvm-cov report ./build-cov/tests/tests \
  -instr-profile=build-cov/tests.profdata \
  -ignore-filename-regex='(extern|tests|build)/'

# 4b. Line-level annotation for one file
xcrun llvm-cov show ./build-cov/tests/tests \
  -instr-profile=build-cov/tests.profdata \
  xdiag/symmetries/tables/representative_table.cpp

# 4c. HTML report
xcrun llvm-cov show ./build-cov/tests/tests \
  -instr-profile=build-cov/tests.profdata \
  -format=html -output-dir=misc/coverage/coverage-html \
  -ignore-filename-regex='(extern|tests|build)/'
```

> On macOS, use the Xcode toolchain's `llvm-cov` / `llvm-profdata` (via `xcrun`,
> as the script does) so their version matches the Apple Clang that produced the
> data. Do **not** mix in a Homebrew `llvm-cov` — version mismatches typically
> fail to read the coverage mapping. On Linux, plain `llvm-cov` / `llvm-profdata`
> that match your `clang` are fine.

## Notes specific to xdiag

- **Templates / headers.** Source-based coverage counts each instantiation, so a
  template exercised for only one type shows the other instantiations as
  uncovered. Read `bits/`, kernel, and table reports with that in mind.
- **OpenMP.** Parallel regions are measured correctly. Set `OMP_NUM_THREADS=1`
  for deterministic region counts while chasing a specific gap.
- **MPI (`tests_distributed`).** Give each rank its own profile via a `%p` in the
  filename, then merge them all:
  ```bash
  LLVM_PROFILE_FILE='build-cov/dist-%p.profraw' mpirun -n 4 ./build-cov/tests/tests_distributed
  xcrun llvm-profdata merge -sparse build-cov/dist-*.profraw -o build-cov/dist.profdata
  ```

## Gitignore

The generated `misc/coverage/coverage-html/` directory and the intermediate
`*.profraw` / `*.profdata` files are gitignored — only `coverage.sh` and this
`README.md` are tracked.
