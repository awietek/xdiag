# SPDX-License-Identifier: Apache-2.0
#
# Read/write the human-readable API model as api/<module>.toml.
#
# Each module file lists its entries grouped under comment banners naming the
# original header. Types are plain strings ("int64_t", "arma::vec", "Op");
# arguments are "type name" or "type name = default". Reading uses tomllib;
# writing is hand-rolled so we control the comments.

import json
import tomllib
from pathlib import Path

API_DIR = Path(__file__).with_name("api")
_FIELD_ORDER = ["kind", "name", "record", "const", "deleted", "returns", "args", "fields"]


def quote(s):
    return json.dumps(s)          # JSON string escaping is valid TOML basic-string


def write(entries):
    by_module = {}
    for e in entries:
        by_module.setdefault(e["module"], []).append(e)

    API_DIR.mkdir(exist_ok=True)
    for stale in API_DIR.glob("*.toml"):
        stale.unlink()

    for module, items in by_module.items():
        out = [f'module = {quote(module)}', ""]
        current_header = None
        for e in items:
            if e["header"] != current_header:
                current_header = e["header"]
                bar = "# " + "=" * 60
                out += ["", bar, f"# {current_header}", bar, ""]
            out.append("[[entries]]")
            for key in _FIELD_ORDER:
                if key not in e:
                    continue
                v = e[key]
                if isinstance(v, bool):
                    out.append(f"{key} = {'true' if v else 'false'}")
                elif isinstance(v, list):
                    out.append(f"{key} = [{', '.join(quote(x) for x in v)}]")
                else:
                    out.append(f"{key} = {quote(v)}")
            out.append("")
        (API_DIR / f"{module}.toml").write_text("\n".join(out).rstrip() + "\n")
    return len(by_module)


def _parse_arg(s):
    default = None
    if " = " in s:
        s, default = s.split(" = ", 1)
    typ, name = s.rsplit(" ", 1)
    return {"type": typ, "name": name, "default": default}


def load():
    """Flat list of entries with module attached and args/fields expanded into
    {type, name, default} dicts."""
    entries = []
    for path in sorted(API_DIR.glob("*.toml")):
        data = tomllib.loads(path.read_text())
        module = data["module"]
        for e in data.get("entries", []):
            e["module"] = module
            e["args"] = [_parse_arg(a) for a in e.get("args", [])]
            e["fields"] = [_parse_arg(f) for f in e.get("fields", [])]
            entries.append(e)
    return entries
