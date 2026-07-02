# SPDX-License-Identifier: Apache-2.0
#
# Sidecar overrides for the Julia wrapper generator.
#
# This is the ONE place where non-derivable wrapper semantics live: what to
# exclude, which C++ types bridge to native Julia types, wrapper-struct field
# names, 1-based<->0-based index shifts, template instantiation combos, and
# routing of pointer-fill routines to hand-written fragments.
#
# Everything here is validated against the AST model at generation time
# (see validate() below): a hard error is raised if any override references a
# symbol / argument / type that no longer exists, so overrides cannot silently
# rot when the C++ API changes.

# ---------------------------------------------------------------------------
# Exclusions
# ---------------------------------------------------------------------------

# Header globs excluded wholesale. Distributed blocks link the MPI build, not
# the serial libxdiagjl, so they are out of scope for the Julia wrapper.
EXCLUDE_HEADERS = [
    "*/blocks/distributed/*",
]

# Records not wrapped as opaque cxx_ types.
EXCLUDE_RECORDS = {
    # Distributed (MPI) — excluded with their iterators.
    "xdiag::SpinhalfDistributed", "xdiag::SpinhalfDistributedIterator",
    "xdiag::tJDistributed", "xdiag::tJDistributedIterator",
    "xdiag::ElectronDistributed", "xdiag::ElectronDistributedIterator",
    # OpSum's internal product structure — construction is covered by the
    # Op / + / * API; revisit if users need explicit monomials.
    "xdiag::Monomial", "xdiag::Term",
    # Internal / not user-facing.
    "xdiag::symmetries::SitePermutation",
    "xdiag::Representation::PermutationIrrep",
    "xdiag::PermutationGroup::iterator",
    "xdiag::io::FileTomlHandler",  # proxy for file["key"]; use read_* free fns
}

# Individual symbols / overloads to skip (fully-qualified). Overloads are
# disambiguated by a tuple of argument category strings where needed.
EXCLUDE_SYMBOLS = {
    # ostream operator<< is replaced by to_string.
    "xdiag::operator<<",
}

# Argument categories that make an overload unwrappable -> that overload is
# dropped (the others for the same name are kept). Used for sugar overloads
# with no clean Julia mapping.
SKIP_OVERLOAD_ARG_CATS = {
    "std_composite",   # std::initializer_list, std::__wrap_iter, chrono, ...
}
# std::istream / std::function overloads handled explicitly below.
SKIP_OVERLOAD_ARG_CPPS = {
    "std::istream",
}

# ---------------------------------------------------------------------------
# Native-type bridges — value types mapped to idiomatic Julia, NOT wrapped.
# ---------------------------------------------------------------------------
# The emitter converts at the boundary instead of exposing an opaque cxx_ type.
#   scalar : Scalar  <-> Number (Float64/ComplexF64)
#   coeff  : Coeff   <-> Union{String,Number}
#   vector : Vector  <-> Vector{Float64|ComplexF64}   (wraps arma)
#   matrix : Matrix  <-> Matrix{Float64|ComplexF64}   (wraps arma)
NATIVE_BRIDGES = {
    "xdiag::Scalar": "scalar",
    "xdiag::Coeff": "coeff",
    "xdiag::Vector": "vector",
    "xdiag::Matrix": "matrix",
}

# ---------------------------------------------------------------------------
# Wrapper-struct field names (Julia `struct T; <field>::cxx_T end`).
# Default is cxx_<snake_case(T)>; only exceptions are listed.
# ---------------------------------------------------------------------------
FIELD_NAMES = {
    # All blocks share `cxx_block` so generic Block dispatch works.
    "xdiag::Spinhalf": "cxx_block",
    "xdiag::tJ": "cxx_block",
    "xdiag::Electron": "cxx_block",
    "xdiag::Boson": "cxx_block",
    "xdiag::Fermion": "cxx_block",
    # snake_case would give the wrong name for these:
    "xdiag::OpSum": "cxx_opsum",
    "xdiag::PermutationGroup": "cxx_group",
    "xdiag::FileToml": "cxx_file",
}

# ---------------------------------------------------------------------------
# 1-based <-> 0-based index shifts (CONFIDENT set).
# ---------------------------------------------------------------------------
# Semantics: a parameter/return that denotes a lattice SITE index is 0-based in
# C++ and 1-based in Julia. The emitter subtracts 1 from flagged args and adds 1
# to flagged returns (element-wise for vector/arma types).
#
# ONEBASED_ARGS: {qualified_symbol: [arg_name, ...]}  (applies across overloads)
ONEBASED_ARGS = {
    "xdiag::Op::Op": ["site", "sites"],
    "xdiag::Permutation::Permutation": ["array", "list"],
    "xdiag::PermutationGroup::PermutationGroup": ["matrix"],
    "xdiag::col": ["n"],
    "xdiag::vector": ["n"],
    "xdiag::vectorC": ["n"],
    "xdiag::fill": ["col"],
}
# ONEBASED_RETURN: {qualified_symbol}  (add 1; element-wise for vectors)
ONEBASED_RETURN = {
    "xdiag::Op::sites",
    "xdiag::Op::operator[]",
    "xdiag::Permutation::array",
    # block index() -> position of a ProductState in the basis (1-based)
    "xdiag::Spinhalf::index",
    "xdiag::tJ::index",
    "xdiag::Electron::index",
    "xdiag::Boson::index",
    "xdiag::Fermion::index",
}

# ---------------------------------------------------------------------------
# AMBIGUOUS shifts — NEED A RULING, deliberately NOT applied yet.
# ---------------------------------------------------------------------------
# The current hand-written wrapper is inconsistent here; encoding either way
# changes behaviour, so these are quarantined until confirmed. Each entry:
#   (symbol, what, current_wrapper_behaviour, semantic_expectation)
AMBIGUOUS_ONEBASED_REVIEW = [
    ("xdiag::Op::operator[]", "arg idx",
     "NOT shifted (idx passed as-is), return +1",
     "idx is a 1-based position -> arg should be idx-1"),
    ("xdiag::ProductState::operator[]", "arg i / return",
     "NOT shifted either way; returns a local-state value",
     "i is a 1-based site -> arg i-1; return is a local state, no shift"),
    ("xdiag::ProductState::setindex!", "arg idx",
     "NOT shifted",
     "idx is a 1-based site -> arg idx-1"),
    ("xdiag::Permutation::operator[]", "arg i / return",
     "not wrapped currently",
     "i is a 1-based site -> arg i-1; return is a site -> +1"),
]

# ---------------------------------------------------------------------------
# Template instantiation combos.
# ---------------------------------------------------------------------------
# Concrete type tuples to instantiate for each template symbol. Keyed by
# fully-qualified name; the emitter produces one wrapped method per combo.
# Note: the new (memory-branch) kernels take `Block const&` (the block variant)
# rather than a block_t template parameter, so blocks need NO instantiation.
IDX_TYPES = ["int64_t", "int32_t"]
COEFF_TYPES = ["double", "complex"]
OP_TYPES = ["Op", "OpSum"]           # Monomial excluded

# The Julia wrapper builds all sparse formats from the CSR two-phase functions
# (mirroring the C++ library, where CSC is a transposed CSR): csr_matrix_nnz ->
# allocate Julia-side -> csr_matrix_fill. CSC/COO are derived Julia-side, so
# their nnz/fill are intentionally NOT wrapped here.
TEMPLATES = {
    # dense matrix builders
    "xdiag::matrix":  {"op_t": OP_TYPES},          # matrix(op_t, Block[, Block])
    "xdiag::matrixC": {"op_t": OP_TYPES},
    # sparse two-phase build (CSR only; special-cased below)
    "xdiag::csr_matrix_nnz":  {"coeff_t": COEFF_TYPES},
    "xdiag::csr_matrix_fill": {"idx_t": IDX_TYPES, "coeff_t": COEFF_TYPES},
    # dense conversion of sparse
    "xdiag::to_dense": {"idx_t": IDX_TYPES, "coeff_t": COEFF_TYPES},
}

# Template symbols whose combos are still TODO (resolved in Phase 3 against the
# compiler): apply() over mat_t, and the sparse-matrix-input variants of the
# lanczos / lobpcg / time-evolution routines. Listed so coverage stays honest.
TEMPLATES_TODO = {
    "xdiag::apply",
    "xdiag::eigs_lanczos", "xdiag::eigvals_lanczos", "xdiag::eigvals_lanczos_inplace",
    "xdiag::eigs_lobpcg",
    "xdiag::eig0", "xdiag::eigval0", "xdiag::eigs", "xdiag::eigvals",
    "xdiag::evolve_lanczos", "xdiag::evolve_lanczos_inplace",
    "xdiag::time_evolve", "xdiag::time_evolve_inplace",
    "xdiag::time_evolve_expokit", "xdiag::time_evolve_expokit_inplace",
    "xdiag::imaginary_time_evolve", "xdiag::imaginary_time_evolve_inplace",
    "xdiag::csr_matrix", "xdiag::csc_matrix", "xdiag::coo_matrix",  # concrete + template overloads
}

# ---------------------------------------------------------------------------
# Specials — hand-written fragments spliced verbatim (pointer-fill routines).
# ---------------------------------------------------------------------------
# The Julia side allocates, then calls these to fill its own memory.
#   dense:  matrix<coeff_t>(OpSum, Block, Block, coeff_t* mat)
#   sparse: *_matrix_nnz -> allocate -> *_matrix_fill(..., idx_t*, idx_t*, coeff_t*)
SPECIALS = {
    "xdiag::matrix": "matrix",           # only the coeff_t* pointer overload
    "xdiag::csr_matrix_nnz": "sparse_csr",   # phase 1: per-row counts
    "xdiag::csr_matrix_fill": "sparse_csr",  # phase 2: fill Julia-owned arrays
    # GPWF / fill std::function callbacks: Julia-function -> std::function bridge
    "xdiag::fill": "fill_function",      # only the coeff_f std::function overloads
}


# ---------------------------------------------------------------------------
# Validation against the AST model.
# ---------------------------------------------------------------------------

def _index_model(model):
    by_qual = {}
    arg_names = {}
    for e in model["entries"]:
        q = e.get("qualified")
        if not q:
            continue
        by_qual.setdefault(q, []).append(e)
        for a in e.get("args", []):
            arg_names.setdefault(q, set()).add(a["name"])
    return by_qual, arg_names


def validate(model):
    """Hard-error if any override references something the model doesn't have."""
    by_qual, arg_names = _index_model(model)
    records = set(model["record_names"])
    errors = []

    def need_symbol(q, where):
        if q not in by_qual:
            errors.append(f"{where}: unknown symbol '{q}'")

    def need_record(q, where):
        if q not in records:
            errors.append(f"{where}: unknown record '{q}'")

    for q in EXCLUDE_RECORDS:
        need_record(q, "EXCLUDE_RECORDS")
    for q in NATIVE_BRIDGES:
        need_record(q, "NATIVE_BRIDGES")
    for q in FIELD_NAMES:
        need_record(q, "FIELD_NAMES")

    for q, names in ONEBASED_ARGS.items():
        need_symbol(q, "ONEBASED_ARGS")
        have = arg_names.get(q, set())
        for n in names:
            if have and n not in have:
                errors.append(
                    f"ONEBASED_ARGS['{q}']: no argument named '{n}' "
                    f"(have: {sorted(have)})")
    for q in ONEBASED_RETURN:
        need_symbol(q, "ONEBASED_RETURN")

    for q in SPECIALS:
        need_symbol(q, "SPECIALS")
    for q, combos in TEMPLATES.items():
        need_symbol(q, "TEMPLATES")
        ents = [e for e in by_qual.get(q, []) if e.get("is_template")]
        tparams = set()
        for e in ents:
            tparams.update(e.get("template_params", []))
        for p in combos:
            if tparams and p not in tparams:
                errors.append(
                    f"TEMPLATES['{q}']: no template parameter '{p}' "
                    f"(have: {sorted(tparams)})")

    return errors


if __name__ == "__main__":
    import sys
    sys.path.insert(0, __import__("os").path.dirname(__file__))
    import model as _model
    m = _model.build_model([])
    errs = validate(m)
    if errs:
        print(f"OVERRIDES VALIDATION FAILED ({len(errs)} issue(s)):", file=sys.stderr)
        for e in errs:
            print("  - " + e, file=sys.stderr)
        sys.exit(1)
    print("overrides.py: validation OK — all references resolve against the model.")
