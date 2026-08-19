# MkDocs hook: render a BibTeX file into a formatted, year-sorted publication
# list. Any page containing the token `<!-- bibliography -->` gets it replaced
# by the rendered list (parsed from docs/assets/papers.bib).
#
# Self-contained: uses a minimal BibTeX parser so no external package
# (bibtexparser / pybtex) needs to be installed.

import os
import re

TOKEN = "<!-- bibliography -->"
BIB_RELPATH = os.path.join("assets", "papers.bib")

# Within each year, papers (co)authored by this name are listed last
# (matched case-insensitively against the `author` field).
SELF_AUTHOR = "wietek"


def _clean(value):
    # Strip BibTeX braces and normalize whitespace / common escapes.
    value = value.replace("{", "").replace("}", "")
    value = value.replace("\\&", "&").replace("--", "–")
    return re.sub(r"\s+", " ", value).strip()


def _parse_fields(text):
    fields = {}
    i, n = 0, len(text)
    while i < n:
        m = re.match(r"\s*(\w+)\s*=\s*", text[i:])
        if not m:
            break
        name = m.group(1).lower()
        i += m.end()
        if i >= n:
            break
        if text[i] == "{":
            depth, j = 1, i + 1
            while j < n and depth > 0:
                depth += (text[j] == "{") - (text[j] == "}")
                j += 1
            value, i = text[i + 1 : j - 1], j
        elif text[i] == '"':
            j = text.find('"', i + 1)
            value, i = text[i + 1 : j], j + 1
        else:
            m2 = re.match(r"([^,\n]+)", text[i:])
            value, i = m2.group(1), i + m2.end()
        fields[name] = _clean(value)
        i += re.match(r"\s*,?\s*", text[i:]).end()
    return fields


def _parse_bib(text):
    entries, i, n = [], 0, len(text)
    while True:
        at = text.find("@", i)
        if at == -1:
            break
        m = re.match(r"@(\w+)\s*\{", text[at:])
        if not m:
            i = at + 1
            continue
        start = at + m.end()
        depth, j = 1, start
        while j < n and depth > 0:
            depth += (text[j] == "{") - (text[j] == "}")
            j += 1
        body = text[start : j - 1]
        i = j
        if m.group(1).lower() in ("comment", "string", "preamble"):
            continue
        comma = body.find(",")
        if comma == -1:
            continue
        entry = {"key": body[:comma].strip()}
        entry.update(_parse_fields(body[comma + 1 :]))
        entries.append(entry)
    return entries


def _format_authors(authors):
    parts = [p.strip() for p in re.split(r"\s+and\s+", authors) if p.strip()]
    if len(parts) > 1:
        return ", ".join(parts[:-1]) + ", and " + parts[-1]
    return parts[0] if parts else ""


def _format_entry(e):
    authors = _format_authors(e.get("author", ""))
    title = e.get("title", "")
    year = e.get("year", "")
    eid = re.sub(r"(?i)^arxiv:", "", e["eprint"]) if e.get("eprint") else None

    venue = e.get("journal") or e.get("booktitle") or ""
    if venue:
        if e.get("volume"):
            venue += f" **{e['volume']}**"
        if e.get("pages"):
            venue += f", {e['pages']}"
    elif eid:
        venue = f"arXiv:{eid}"

    url = None
    if e.get("doi"):
        url = f"https://doi.org/{e['doi']}"
    elif e.get("url"):
        url = e["url"]
    elif eid:
        url = f"https://arxiv.org/abs/{eid}"

    line = f"- {authors}, *{title}*" if authors else f"- *{title}*"
    if venue:
        line += f", {venue}"
    if year:
        line += f" ({year})"
    line += "."
    if url:
        line += f" [link]({url})"
    return line


_MONTHS = {
    m: i
    for i, m in enumerate(
        ["jan", "feb", "mar", "apr", "may", "jun",
         "jul", "aug", "sep", "oct", "nov", "dec"], 1)
}


def _year(e):
    try:
        return int(re.sub(r"\D", "", e.get("year", "0")) or 0)
    except ValueError:
        return 0


def _month(e):
    return _MONTHS.get(e.get("month", "").strip().lower()[:3], 0)


def _render(bib_path):
    if not os.path.exists(bib_path):
        return "_No publications listed yet._"
    entries = _parse_bib(open(bib_path, encoding="utf-8").read())
    if not entries:
        return "_No publications listed yet._"

    # Sort: newest year first; within a year, papers not (co)authored by
    # SELF_AUTHOR come first and self-authored ones last; then by month.
    entries.sort(
        key=lambda e: (
            -_year(e),
            SELF_AUTHOR in e.get("author", "").lower(),
            -_month(e),
        )
    )

    # group under a heading per year
    blocks, current = [], None
    for e in entries:
        year = e.get("year", "").strip() or "Unknown"
        if year != current:
            current = year
            blocks.append(f"\n## {year}\n")
        blocks.append(_format_entry(e))
    return "\n".join(blocks).strip()


def on_page_markdown(markdown, page, config, files, **kwargs):
    if TOKEN not in markdown:
        return markdown
    bib_path = os.path.join(config["docs_dir"], BIB_RELPATH)
    return markdown.replace(TOKEN, _render(bib_path))
