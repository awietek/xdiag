#!/bin/bash
# SPDX-License-Identifier: Apache-2.0
#
# Bootstrap a libclang venv (first run only) and regenerate the Markdown
# reference of the public xdiag C++ API into misc/api.md.
#
# Usage:  misc/generate_api_md.sh [extra clang args...]

set -e
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(dirname "$script_dir")"
venv="$script_dir/.venv"

if [ ! -x "$venv/bin/python" ]; then
    echo "Setting up libclang venv in $venv ..."
    python3 -m venv "$venv"
    "$venv/bin/pip" install -q --upgrade pip libclang
fi

"$venv/bin/python" "$script_dir/generate_api_md.py" \
    -o "$script_dir/api.md" "$@"
