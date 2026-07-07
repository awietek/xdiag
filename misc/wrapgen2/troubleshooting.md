# XDiag.jl: "Unregistered method with key" on `using XDiag`

## Symptom

```
ERROR: InitError: Unregistered method with key
(:constructor, :cxx_arma_vec, :XDiag, 0x900757a253af782e) requested,
maybe you need to precompile the Julia module?
```

Thrown from `@initcxx` → `CxxWrap.initialize_julia_module` during `using XDiag`.

## Cause

A **stale Julia precompile cache** — not a broken library.

- `@wrapmodule` reads the C++ library **at precompile time** and bakes the C++
  method keys into the `.ji` cache.
- `@initcxx` re-binds those keys against the **live** library at runtime.
- Julia's cache-staleness check tracks only Julia source files and dependency
  *package versions*. It does **not** watch the content of an external
  `.dylib`.
- `XDiag_jll` is overridden (`~/.julia/artifacts/Overrides.toml`) to point at
  `.../xdiag/install`. So rebuilding `libxdiagjl.dylib` swaps the library
  underneath Julia **without touching any `.jl` source** → the old `.ji` still
  looks valid and gets deserialized.
- The old cache carried the previous ABI name `cxx_arma_vec` (from when
  `XDiag.jl` still did `include("utils/armadillo.jl")`); the rebuilt lib
  registers the new wrapgen2 name `typ_arma_vec`. Key mismatch → error.

## Fix

After every rebuild + reinstall of `libxdiagjl`:

```bash
rm -rf ~/.julia/compiled/*/XDiag
```

Then `using XDiag` recompiles cleanly against the current library.

## Make it routine

Add to your build/install script so it can never drift:

```bash
cmake --build build --target install   # rebuild + reinstall libxdiagjl
rm -rf ~/.julia/compiled/*/XDiag        # invalidate stale precompile caches
```

## Diagnostic one-liners

```bash
# Which ABI name is in the installed lib?
strings install/lib/libxdiagjl.dylib | grep -o 'cxx_arma_vec\|typ_arma_vec' | sort | uniq -c

# Which ABI name is baked into each cached .ji?
for f in ~/.julia/compiled/*/XDiag/*.ji; do
  echo "$f"; strings "$f" | grep -o 'cxx_arma_vec\|typ_arma_vec' | sort | uniq -c
done

# Where does JLL resolve? (active line = uncommented)
cat ~/.julia/artifacts/Overrides.toml

# Which XDiag source does the current Julia env use?
julia -e 'println(Base.find_package("XDiag"))'
```

## Environment notes / traps

- Default `julia` (1.12.x) `dev`s the new minimal wrapper at
  `~/.julia/dev/XDiag`. The v1.11 env resolves to the **old registry copies**
  in `~/.julia/packages/XDiag/{uKdat,DcSj0}` (full package, old
  `cxx_arma_vec` ABI) — don't test the new lib there.
- **Latent mismatch:** `~/.julia/dev/XDiag/src/utils/armadillo.jl` still uses
  the old names (`cxx_arma_vec`, `cxx_arma_mat`, `cxx_arma_ivec`), but the new
  lib registers `typ_arma_vec` / `typ_arma_mat` /
  `typ_arma_vec_int64_t`. Reconcile these (rename in the C++ `add_type<>` calls
  or add Julia-side aliases) before re-enabling the ergonomic layer in
  `XDiag.jl`.
