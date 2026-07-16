# SPDX-License-Identifier: Apache-2.0
#
# The only non-derivable policy: 1-based<->0-based index shifts, and the
# handful of symbols we do NOT generate (provided by hand-written static files,
# or genuinely non-mapping). Everything is keyed by (record, name); free
# functions use record "". validate() errors if a key no longer exists.

# args whose value is a site/element index: -1 into C++, +1 back out
ONEBASED_ARGS = {
    ("Op", "Op"): ["site", "sites"],
    ("Op", "operator[]"): ["idx"],
    ("Permutation", "Permutation"): ["array"],
    ("Permutation", "operator[]"): ["i"],
    ("ProductState", "operator[]"): ["i"],
    ("PermutationGroup", "operator[]"): ["sym"],
    # State column accessors: n selects the n-th vector (1-based in Julia)
    ("", "vector"): ["n"],
    ("", "vectorC"): ["n"],
    ("", "col"): ["n"],
    ("State", "col"): ["n"],
}
# returns that are themselves a site/element index: +1
ONEBASED_RETURN = {
    ("Op", "sites"), ("Op", "operator[]"),
    ("Permutation", "array"), ("Permutation", "operator[]"),
    ("Spinhalf", "index"), ("tJ", "index"), ("Electron", "index"),
    ("Boson", "index"), ("Fermion", "index"),
}
# Whole classes never wrapped (not add_type'd; args/returns using them skip).
IGNORED_CLASSES = {"FileTomlHandler"}

# Deliberately out of scope: never wrapped, never hand-written.
IGNORED = {
    "begin", "end", "what",
    "operator<<", "operator<", "operator<=", "operator>", "operator>=",
    "operator!=", "operator++",
    "empty",   # C++-internal only; collides with Base.empty and unneeded in Julia
}
# Excluded from generation because a hand-written static file provides them
# (e.g. the pointer-fill matrix builders).
HAND_WRITTEN = {"fill"}

# Iteration via the C++ iterator protocol (begin()/end() + iterator_t with
# ++, unary *, ==). container -> element wrapper. Blocks yield ProductStates.
ITERABLES = {
    "Spinhalf": "ProductState", "tJ": "ProductState", "Electron": "ProductState",
    "Boson": "ProductState", "Fermion": "ProductState",
}
# Iteration by size() + getindex (no extra C++); the element is whatever
# getindex returns. ProductState yields its local states as Int64.
INDEXED_ITERABLES = {"ProductState"}


def key(e):
    return (e.get("record") or "", e["name"])


def arg_shifts(e):
    return set(ONEBASED_ARGS.get(key(e), []))


def return_shift(e):
    return key(e) in ONEBASED_RETURN


def skip_reason(e):
    """'ignored' / 'hand-written' / None -- why a symbol is not generated."""
    if e["name"] in IGNORED:
        return "ignored"
    if e["name"] in HAND_WRITTEN:
        return "hand-written"
    return None


def validate(entries):
    symbols = {key(e) for e in entries}
    arg_names = {}
    for e in entries:
        arg_names.setdefault(key(e), set()).update(a["name"] for a in e.get("args", []))
    problems = []
    for k, names in ONEBASED_ARGS.items():
        if k not in symbols:
            problems.append(f"ONEBASED_ARGS: no such symbol {k}")
        else:
            for n in names:
                if n not in arg_names.get(k, set()):
                    problems.append(f"ONEBASED_ARGS{k}: no argument '{n}'")
    for k in ONEBASED_RETURN:
        if k not in symbols:
            problems.append(f"ONEBASED_RETURN: no such symbol {k}")
    return problems
