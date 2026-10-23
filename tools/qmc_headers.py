"""Minimal C header parser used by tools/gen_api.py

Not a general C parser: it understands exactly declaration shapes used in
QMC's public headers (prototypes, typedef struct/enum/union, function-pointer
typedefs, static inline functions, #define macros, extern data) and attaches
comment that precedes each one
"""
import re

IDENT = re.compile(r"[A-Za-z_]\w*")


class Decl:
    def __init__(self, kind, name, text, comment, line):
        self.kind = kind          # 'function' | 'type' | 'macro' | 'data'
        self.name = name
        self.text = text          # source text of declaration
        self.comment = comment    # cleaned preceding comment ('' if none)
        self.line = line


class Header:
    def __init__(self, path):
        self.path = path
        self.module_doc = ""      # file-level description ('' if none)
        self.decls = []


def clean_comment(raw):
    raw = raw.strip()
    if raw.startswith("/*"):
        raw = raw[2:]
        if raw.endswith("*/"):
            raw = raw[:-2]
        lines = []
        for ln in raw.split("\n"):
            ln = ln.strip()
            ln = re.sub(r"^\*+ ?", "", ln)
            lines.append(ln.rstrip())
    else:  # run of // lines
        lines = [re.sub(r"^//+ ?", "", ln.strip()) for ln in raw.split("\n")]
    while lines and not lines[0].strip():
        lines.pop(0)
    while lines and not lines[-1].strip():
        lines.pop()
    return "\n".join(lines)


def tokenize(src):
    """Yield ('comment'|'pp'|'code', text, start_line, end_line, blank_after)."""
    i, n, line = 0, len(src), 1
    out = []
    code_start = None
    code_line = 0

    def flush_code(end):
        nonlocal code_start
        if code_start is not None:
            txt = src[code_start:end]
            if txt.strip():
                out.append(["code", txt.strip(), code_line, line, False])
            code_start = None

    while i < n:
        c = src[i]
        if c == "\n":
            line += 1
            i += 1
            continue
        if src.startswith("/*", i) and (code_start is None):
            j = src.find("*/", i + 2)
            j = n if j < 0 else j + 2
            start_line = line
            line += src.count("\n", i, j)
            out.append(["comment", src[i:j], start_line, line, False])
            i = j
            continue
        if src.startswith("//", i) and (code_start is None):
            j = i
            start_line = line
            buf = []
            while True:
                k = src.find("\n", j)
                k = n if k < 0 else k
                buf.append(src[j:k])
                j = k + 1
                line += 1
                m = re.match(r"[ \t]*//", src[j:j + 200]) if j < n else None
                if not m:
                    break
                j += m.end() - 2
            out.append(["comment", "\n".join(buf), start_line, line - 1, False])
            i = min(j, n)
            continue
        if c == "#" and (code_start is None) and (i == 0 or src[i - 1] == "\n"):
            j = i
            start_line = line
            while True:
                k = src.find("\n", j)
                k = n if k < 0 else k
                seg = src[j:k]
                j = k + 1
                line += 1
                if not seg.rstrip().endswith("\\"):
                    break
            out.append(["pp", src[i:j].rstrip(), start_line, line - 1, False])
            i = min(j, n)
            continue
        if c in " \t\r":
            i += 1
            continue
        m = re.match(r'extern\s+"C"\s*\{', src[i:i + 40])
        if m:                      # `extern "C" {` opens no statement
            i += m.end()
            continue
        if c == "}":               # its closing brace
            i += 1
            continue
        # code: scan to statement end, respecting braces/parens and comments
        code_start = i
        code_line = line
        depth = 0
        j = i
        while j < n:
            ch = src[j]
            if ch == "\n":
                line += 1
            if src.startswith("/*", j):
                k = src.find("*/", j + 2)
                k = n if k < 0 else k + 2
                line += src.count("\n", j, k)
                j = k
                continue
            if src.startswith("//", j):
                k = src.find("\n", j)
                j = n if k < 0 else k
                continue
            if ch in "({[":
                depth += 1
            elif ch in ")}]":
                depth -= 1
                if ch == "}" and depth == 0:
                    # function body or struct body ends; a struct decl continues
                    # to its ';' ; a function body ends statement
                    head = src[code_start:j]
                    rest = src[j + 1:j + 200]
                    if not re.match(r"\s*[A-Za-z_]\w*\s*;|\s*;|\s*\*?\s*[A-Za-z_]", rest) or \
                            re.match(r"\s*(?:static|inline|typedef)", head) is None and ")" in head.split("{")[0][-3:] + ")":
                        pass
            elif ch == ";" and depth == 0:
                j += 1
                break
            if ch == "}" and depth == 0:
                # statement ending with a closing brace and no ';' => function body
                k = j + 1
                tail = src[k:k + 80]
                if not re.match(r"\s*(?:[A-Za-z_]\w*\s*)?(?:,|;)", tail):
                    j += 1
                    break
            j += 1
        out.append(["code", src[code_start:j].strip(), code_line, line, False])
        code_start = None
        i = j
    return out


def parse_header(path):
    src = open(path, encoding="utf-8").read()  # noqa: SIM115
    h = Header(path)
    toks = tokenize(src)
    pending = None
    pending_end = -10
    seen_code = False
    for idx, (kind, text, l0, l1, _) in enumerate(toks):
        nxt = toks[idx + 1] if idx + 1 < len(toks) else None
        if kind == "comment":
            pending = clean_comment(text)
            pending_end = l1
            # a standalone (blank line after) block before first real
            # declaration is module description
            if not seen_code and nxt is not None and \
                    (nxt[2] - l1 > 1 or nxt[0] == "pp") and \
                    len(pending.split()) >= 8 and not h.module_doc:
                h.module_doc = pending
                pending = None
            continue
        if kind == "pp":
            m = re.match(r"#\s*define\s+(\w+)(\(.*?\))?\s*(.*)", text, re.S)
            if m and not re.match(r"QMC_\w+_H$", m.group(1)) and \
                    not m.group(1).endswith("_H"):
                comment = pending if pending and l0 - pending_end <= 1 else ""
                h.decls.append(Decl("macro", m.group(1), text, comment, l0))
            if not re.match(r"#\s*(ifndef|define|endif|if|ifdef|else|elif|include|pragma|undef)", text):
                pass
            if not text.lstrip("#").strip().startswith(("ifndef", "define")):
                pass
            if text.lstrip("# ").startswith("include"):
                pass
            pending = pending if (pending and l0 - pending_end <= 1) else None
            continue
        # code
        if text.startswith('extern "C"') or text in ("}", "} /* extern C */"):
            continue
        t = text.strip()
        if not t or t == "}":
            continue
        comment = pending if pending and l0 - pending_end <= 2 else ""
        if not h.module_doc and not seen_code and comment and \
                len(comment.split()) >= 12 and pending is not None and False:
            h.module_doc = comment
        seen_code = True
        pending = None
        d = classify(t, comment, l0)
        if d:
            h.decls.append(d)
    return h


def strip_comments(s):
    s = re.sub(r"/\*.*?\*/", " ", s, flags=re.S)
    s = re.sub(r"//[^\n]*", " ", s)
    return s


def classify(t, comment, line):
    bare = strip_comments(t)
    b = bare.strip()
    if b.startswith("typedef"):
        # typedef struct {...} name;  /  typedef enum {...} name;  /  fn ptr
        m = re.search(r"\(\s*\*\s*(\w+)\s*\)\s*\(", b)
        if m:
            return Decl("type", m.group(1), t, comment, line)
        m = re.search(r"(\w+)\s*;\s*$", b)
        return Decl("type", m.group(1) if m else "?", t, comment, line)
    if re.match(r"(struct|enum|union)\b", b):
        m = re.match(r"(?:struct|enum|union)\s+(\w+)", b)
        return Decl("type", m.group(1) if m else "?", t, comment, line)
    if "(" in b:
        head = b.split("{")[0]
        m = re.search(r"([A-Za-z_]\w*)\s*\(", head)
        # skip pure function-pointer variable decls, attribute macros
        if m:
            name = m.group(1)
            # handle QMC_INTERNAL-style macro prefixes: take last ident before '('
            idents = re.findall(r"([A-Za-z_]\w*)\s*\(", head)
            name = idents[0]
            sig = head.strip().rstrip(";").strip()
            return Decl("function", name, sig + ";", comment, line)
    m = re.search(r"(\w+)\s*(?:\[[^\]]*\])?\s*(?:=[^;]*)?;\s*$", b)
    if m:
        return Decl("data", m.group(1), t, comment, line)
    return None
