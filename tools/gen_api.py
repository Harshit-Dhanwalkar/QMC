#!/usr/bin/env python3
"""Generate QMC's API reference, module catalog and umbrella-header block

  python3 tools/gen_api.py           regenerate everything
  python3 tools/gen_api.py --check   exit 1 if any generated file is stale

Inputs : public headers (core/, physics/, export/, latex/), one-line
         summaries/categories in tools/api_meta.py and existing
         examples/, tests/ and docs/src/ trees
Outputs: docs/src/api/*.md          per-module API reference (mdBook)
         docs/src/SUMMARY.md        "API Reference" section (everything
                                    from its heading to end of file)
         core/catalog_data.inc      table behind qmc_help()/qmc_module_find()
         qmc.h                      block between gen_api markers

Run it after adding or changing a public header; CI runs --check
"""
import glob
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tools"))
from qmc_headers import parse_header  # noqa: E402
import api_meta  # noqa: E402

# Headers deliberately left out of <qmc/qmc.h>: they need optional third-party
# libraries or a plotting back end.
NOT_IN_UMBRELLA = {"hdf5_writer", "netcdf_writer", "plot"}
HEADER_GLOBS = ["core/*.h", "core/*/*.h", "physics/*.h", "export/*.h", "latex/*.h"]
SKIP_SUFFIX = ("_internal.h",)

CATEGORY_ORDER = [
    api_meta.SP, api_meta.ATOMS, api_meta.OPEN, api_meta.CHEM, api_meta.QMCARLO,
    api_meta.CM, api_meta.QC, api_meta.EXPORT, api_meta.CORE_TYPES,
    api_meta.CORE_LA, api_meta.CORE_FFT, api_meta.CORE_ODE, api_meta.CORE_SPECIAL,
]
MIN_GUIDE_WORDS = 40


def rel(*p):
    return os.path.join(ROOT, *p)


def read(path):
    with open(path, encoding="utf-8") as f:
        return f.read()


def discover():
    heads = []
    for g in HEADER_GLOBS:
        for p in sorted(glob.glob(rel(g))):
            r = os.path.relpath(p, ROOT)
            if r.endswith(SKIP_SUFFIX) or "/gr/" in r:
                continue
            heads.append(r)
    return heads


def stem_of(path):
    return os.path.basename(path)[:-2]


def first_sentence(text):
    text = " ".join(text.split())
    m = re.search(r"(.+?[.!?])(\s|$)", text)
    return (m.group(1) if m else text)


def find_guide(stem):
    cands = []
    if stem in api_meta.GUIDES and api_meta.GUIDES[stem]:
        cands.append(api_meta.GUIDES[stem])
    cands += ["physics/%s.md" % stem, "internals/%s.md" % stem]
    for c in cands:
        p = rel("docs", "src", c)
        if os.path.exists(p) and len(read(p).split()) >= MIN_GUIDE_WORDS:
            return c
    return ""


def scan_usage():
    """Map header stem -> (examples, tests) that include it or are named after it."""
    ex = sorted(glob.glob(rel("examples", "*.c")))
    te = sorted(glob.glob(rel("tests", "test_*.c")))
    inc = re.compile(r'#\s*include\s+[<"](?:[^">]*/)?([A-Za-z0-9_]+)\.h[">]')
    included = {}
    for p in ex + te:
        r = os.path.relpath(p, ROOT)
        for m in inc.finditer(read(p)):
            included.setdefault(m.group(1), []).append(r)
    return ex, te, included


def usage_for(stem, ex, te, included, is_core):
    def named(paths, pat):
        return [os.path.relpath(p, ROOT) for p in paths
                if re.search(pat, os.path.basename(p))]
    examples = named(ex, r"^eg_\d+_%s\.c$" % re.escape(stem)) or \
        named(ex, r"%s" % re.escape(stem))
    tests = named(te, r"^test_%s\.c$" % re.escape(stem)) or \
        named(te, r"%s" % re.escape(stem))
    if not is_core:
        users = included.get(stem, [])
        if not examples:
            examples = [u for u in users if u.startswith("examples/")][:3]
        if not tests:
            tests = [u for u in users if u.startswith("tests/")][:3]
    return examples[:3], tests[:3]


def sig_oneline(sig):
    return " ".join(sig.split())


def build_modules():
    heads = discover()
    ex, te, included = scan_usage()
    mods = []
    missing = []
    seen = {}
    for h in heads:
        stem = stem_of(h)
        if stem in seen:
            raise SystemExit("duplicate module stem %r: %s and %s" % (stem, seen[stem], h))
        seen[stem] = h
        if stem not in api_meta.SUMMARIES:
            missing.append(h)
            continue
        cat, summary = api_meta.SUMMARIES[stem]
        if len(summary) > 140:
            raise SystemExit("summary for %s is %d chars (max 140)" % (stem, len(summary)))
        parsed = parse_header(rel(h))
        funcs = [d.name for d in parsed.decls if d.kind == "function"]
        start = api_meta.START.get(stem)
        if start:
            for s in start:
                if s not in funcs:
                    raise SystemExit("api_meta.START[%r] names unknown function %s" % (stem, s))
        else:
            # entry points: skip lifecycle/accessor helpers (free, close, ...)
            boring = re.compile(r"(_free|_destroy|_close|_strerror|_default|_dim|_orbitals|_count)$")
            start = [f for f in funcs if not boring.search(f)][:5] or funcs[:5]
        is_core = h.startswith("core/")
        examples, tests = usage_for(stem, ex, te, included, is_core)
        mods.append(dict(
            stem=stem, header=h, category=cat, summary=summary, parsed=parsed,
            start=start, guide=find_guide(stem),
            example=examples[0] if examples else "", examples=examples,
            test=tests[0] if tests else "", tests=tests,
            in_umbrella=0 if stem in NOT_IN_UMBRELLA else 1,
        ))
    if missing:
        raise SystemExit(
            "no entry in tools/api_meta.py SUMMARIES for: " + ", ".join(missing) +
            "\nAdd (category, one-line summary) for each, then rerun.")
    stale = sorted(set(api_meta.SUMMARIES) - set(seen))
    if stale:
        raise SystemExit("api_meta.SUMMARIES lists modules with no header: " + ", ".join(stale))
    order = {c: i for i, c in enumerate(CATEGORY_ORDER)}
    for m in mods:
        if m["category"] not in order:
            raise SystemExit("unknown category %r for %s" % (m["category"], m["stem"]))
    mods.sort(key=lambda m: (order[m["category"]], m["stem"]))
    return mods


# Pages
def md_escape_cell(s):
    return s.replace("|", "\\|")


def render_module(m, all_mods):
    p = m["parsed"]
    stem, header = m["stem"], m["header"]
    funcs = [d for d in p.decls if d.kind == "function"]
    types = [d for d in p.decls if d.kind == "type"]
    macros = [d for d in p.decls if d.kind == "macro"]
    data = [d for d in p.decls if d.kind == "data"]
    out = []
    out.append("# `%s`\n" % stem)
    out.append("> Generated from `%s` by `tools/gen_api.py`. To change this page, "
               "edit header's comments or `tools/api_meta.py`.\n" % header)
    out.append("%s\n" % m["summary"])
    out.append("**Category:** %s\n" % m["category"])
    out.append("```c\n#include <qmc/%s>%s\n```\n" % (
        header, "" if m["in_umbrella"] else
        "   /* not part of <qmc/qmc.h>: include it directly */"))
    if m["start"]:
        out.append("**Start with:** " + ", ".join("[`%s`](#%s)" % (n, n.lower())
                                                    for n in m["start"]) + "\n")
    rows = []
    if m["guide"]:
        rows.append("- Guide: [%s](../%s)" % (m["guide"], m["guide"]))
    for e in m["examples"]:
        rows.append("- Example: `%s`" % e)
    for t in m["tests"]:
        rows.append("- Test (every call with expected results): `%s`" % t)
    rows.append("- Runtime help: `qmc_help(stdout, \"%s\");` or `qmc-help %s`" % (stem, stem))
    out.append("\n".join(rows) + "\n")
    if p.module_doc:
        out.append("## Overview\n")
        out.append("```text\n%s\n```\n" % p.module_doc)
    if types:
        out.append("## Types\n")
        for d in types:
            out.append("### `%s`\n" % d.name)
            if d.comment:
                out.append("```text\n%s\n```\n" % d.comment)
            out.append("```c\n%s\n```\n" % d.text.strip())
    if macros:
        out.append("## Macros and constants\n")
        for d in macros:
            line = d.text.strip()
            out.append("```c\n%s\n```\n" % line)
            if d.comment:
                out.append("```text\n%s\n```\n" % d.comment)
    if data:
        out.append("## Data\n")
        for d in data:
            out.append("```c\n%s\n```\n" % d.text.strip())
    if funcs:
        out.append("## Functions (%d)\n" % len(funcs))
        out.append(", ".join("[`%s`](#%s)" % (d.name, d.name.lower()) for d in funcs) + "\n")
        for d in funcs:
            out.append("### `%s`\n" % d.name)
            out.append("```c\n%s\n```\n" % d.text.strip())
            if d.comment:
                out.append("```text\n%s\n```\n" % d.comment)
            else:
                out.append("_No description in header yet._\n")
    return "\n".join(out).rstrip() + "\n"


def render_index(mods):
    nf = sum(1 for m in mods for d in m["parsed"].decls if d.kind == "function")
    out = ["# API Reference\n",
           "> Generated by `tools/gen_api.py` from the public headers. Do not edit.\n",
           "QMC has **%d modules** and **%d public functions**. Three ways to find "
           "what you need:\n" % (len(mods), nf),
           "1. **This reference**: every module below lists its types, constants and "
           "functions with the comment from the header, plus the example and test that "
           "show it in use.",
           "2. **In code**: `#include <qmc/qmc.h>`, then `qmc_help(stdout, NULL);` for an "
           "overview, `qmc_help(stdout, \"tight_binding\");` for one module, or "
           "`qmc_help(stdout, \"chern\");` to search. "
           "`qmc_module_find()`, `qmc_module_search()` and `qmc_module_at()` give the same "
           "catalog as data.",
           "3. **In a shell**: `qmc-help <module | category | keyword>` (installed by "
           "`make install`).\n",
           "Most functions carry a module prefix (`tb_`, `qaoa_`, `json_`, ...). "
           "See [Using QMC as a Library](../setup/library.md) for installing and linking.\n"]
    cur = None
    for m in mods:
        if m["category"] != cur:
            cur = m["category"]
            out.append("\n## %s\n" % cur)
            out.append("| Module | What it does | Header |")
            out.append("|--------|--------------|--------|")
        out.append("| [`%s`](%s.md) | %s | `%s` |" % (
            m["stem"], m["stem"], md_escape_cell(m["summary"]), m["header"]))
    return "\n".join(out).rstrip() + "\n"


def render_summary(old, mods):
    marker = "# API Reference"
    i = old.find("\n" + marker)
    base = old[:i + 1] if i >= 0 else old.rstrip("\n") + "\n\n"
    out = [base.rstrip("\n") + "\n\n---\n\n" + marker + "\n\n- [Overview](api/index.md)"]
    cur = None
    for m in mods:
        if m["category"] != cur:
            cur = m["category"]
            out.append("- [%s]()" % cur)
        out.append("  - [%s](api/%s.md)" % (m["stem"], m["stem"]))
    text = "\n".join(out) + "\n"
    # collapse a doubled separator if the previous section already ended with one
    return re.sub(r"---\n\n---\n", "---\n", text)


# Catalog
def c_str(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def render_catalog(mods):
    cats = [c for c in CATEGORY_ORDER if any(m["category"] == c for m in mods)]
    out = ["/* GENERATED by tools/gen_api.py from the public headers and tools/api_meta.py.",
           " * Do not edit; run `make api-docs`. Included by core/catalog.c only. */",
           "",
           "static const char *const qmc_catalog_categories[] = {"]
    out += ["    %s," % c_str(c) for c in cats]
    out += ["};", "",
            "static const qmc_module_t qmc_catalog_modules[] = {"]
    for m in mods:
        out.append("    {%s, %s, %s,\n     %s,\n     %s,\n     %s, %s, %s, %d}," % (
            c_str(m["stem"]), c_str(m["header"]), c_str(m["category"]),
            c_str(m["summary"]), c_str(", ".join(m["start"])),
            c_str(m["guide"]), c_str(m["example"]), c_str(m["test"]), m["in_umbrella"]))
    out += ["};", ""]
    return "\n".join(out)


# qmc.h
BEGIN = "/* gen_api:begin (generated by tools/gen_api.py; do not edit) */"
END = "/* gen_api:end */"


def render_umbrella_block(mods):
    lines = [BEGIN]
    groups = [("core", "core/"), ("physics", "physics/"), ("export", "export/"), ("latex", "latex/")]
    for _, pre in groups:
        for h in sorted(m["header"] for m in mods
                        if m["header"].startswith(pre) and m["in_umbrella"] and
                        m["stem"] not in ("fft_wrapper",)):
            lines.append('#include "%s"' % h)
    lines.append(END)
    return "\n".join(lines)


def render_qmc_h(old, mods):
    a, b = old.find(BEGIN), old.find(END)
    if a < 0 or b < 0:
        raise SystemExit("qmc.h is missing the gen_api begin/end markers")
    return old[:a] + render_umbrella_block(mods) + old[b + len(END):]


# Main
def main():
    check = "--check" in sys.argv
    mods = build_modules()
    files = {}
    files["docs/src/api/index.md"] = render_index(mods)
    for m in mods:
        files["docs/src/api/%s.md" % m["stem"]] = render_module(m, mods)
    files["docs/src/SUMMARY.md"] = render_summary(read(rel("docs/src/SUMMARY.md")), mods)
    files["core/catalog_data.inc"] = render_catalog(mods)
    files["qmc.h"] = render_qmc_h(read(rel("qmc.h")), mods)

    stale = []
    for path, text in files.items():
        full = rel(path)
        cur = read(full) if os.path.exists(full) else None
        if cur != text:
            stale.append(path)
            if not check:
                os.makedirs(os.path.dirname(full), exist_ok=True)
                with open(full, "w", encoding="utf-8") as f:
                    f.write(text)
    # pages for modules that no longer exist
    keep = {os.path.basename(p) for p in files if p.startswith("docs/src/api/")}
    for p in glob.glob(rel("docs/src/api/*.md")):
        if os.path.basename(p) not in keep:
            stale.append(os.path.relpath(p, ROOT) + " (removed)")
            if not check:
                os.remove(p)

    if check:
        if stale:
            print("generated API files are out of date:")
            for s in stale:
                print("  " + s)
            print("run `make api-docs` (python3 tools/gen_api.py) and commit result")
            return 1
        print("API reference, catalog and qmc.h are up to date (%d modules)" % len(mods))
        return 0
    print("%d modules; %d file(s) updated" % (len(mods), len(stale)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
