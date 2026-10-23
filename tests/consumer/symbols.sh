#!/bin/sh
# Audit the installed shared library's public surface
#   usage: sh tests/consumer/symbols.sh <prefix>
#
# Fails when libqmc.so
#   1. has no SONAME of the form libqmc.so.<MAJOR> matching core/qmc_version.h
#   2. depends on a sanitizer/coverage runtime (a dev build was installed)
#   3. exports a symbol that is declared in no installed header
#      Such a function is accidentally public: make it `static`, or declare it
#      (Symbols that are intentionally exported but undeclared are listed in
#      tests/consumer/undeclared_symbols_allowed.txt)
set -eu

PREFIX=${1:?usage: symbols.sh <prefix>}
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
LIB="$PREFIX/lib/libqmc.so"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
fail=0

[ -e "$LIB" ] || { echo "FAIL: $LIB not installed"; exit 1; }

major=$(sed -n 's/^#define QMC_VERSION_MAJOR \([0-9]*\)/\1/p' "$ROOT/core/qmc_version.h")
echo "== SONAME"
if command -v readelf >/dev/null 2>&1; then
  soname=$(readelf -d "$LIB" | sed -n 's/.*SONAME.*\[\(.*\)\]/\1/p')
  if [ "$soname" != "libqmc.so.$major" ]; then
    echo "FAIL: SONAME is '$soname', expected libqmc.so.$major"; fail=1
  fi
  echo "== runtime dependencies"
  if readelf -d "$LIB" | grep NEEDED | grep -E 'asan|ubsan|tsan|gcov' >/dev/null; then
    echo "FAIL: libqmc.so needs a sanitizer/coverage runtime:"
    readelf -d "$LIB" | grep NEEDED | grep -E 'asan|ubsan|tsan|gcov'; fail=1
  fi
else
  echo "(readelf not found; skipping SONAME and dependency checks)"
fi

echo "== exported symbols are declared in installed headers"
# Identifier set of every installed header, comments stripped
find "$PREFIX/include/qmc" -name '*.h' | while read -r h; do
  sed -e 's://.*$::' -e ':a;s:/\*[^*]*\*/::g;ta' "$h"
done | tr -cs 'A-Za-z0-9_' '\n' | sort -u > "$TMP/declared.txt"

# Defined, global, non-underscore-prefixed symbols (compiler/linker internals
# such as _init, __bss_start, and OpenMP reduction variables like
# .gomp_critical_user_.reduction.var start with '_' or '.')
nm -D --defined-only "$LIB" | awk '$2 ~ /^[TDBRVW]$/ {print $3}' \
  | grep -v -E '^(_|\.)' | sort -u > "$TMP/exported.txt"

# Intentionally undeclared exports. grep exits 1 when list holds only comments,
# which is normal, empty case
ALLOWED="$HERE/undeclared_symbols_allowed.txt"
[ -f "$ALLOWED" ] || { echo "FAIL: $ALLOWED is missing"; exit 1; }
{ grep -v '^#' "$ALLOWED" || true; } | sort -u > "$TMP/allowed.txt"

undeclared=$(comm -23 "$TMP/exported.txt" "$TMP/declared.txt" | comm -23 - "$TMP/allowed.txt")
if [ -n "$undeclared" ]; then
  echo "FAIL: exported but not declared in any installed header:"
  echo "$undeclared" | sed 's/^/  /'; fail=1
fi
echo "   $(wc -l < "$TMP/exported.txt" | tr -d ' ') exported symbols"

[ "$fail" -eq 0 ] || exit 1
echo "symbol audit passed"
