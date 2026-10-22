#!/bin/sh
# Verify an installed QMC from a downstream project's point of view
#   usage: sh tests/consumer/run.sh <prefix> [compiler]
# Checks:
#   1. every installed public header compiles standalone via <qmc/...>
#   2. pkg-config shared link + run, (3) static link + run
set -eu

PREFIX=${1:?usage: run.sh <prefix> [compiler]}
CC=${2:-cc}
HERE=$(cd "$(dirname "$0")" && pwd)
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

PKG_CONFIG_PATH="$PREFIX/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
export PKG_CONFIG_PATH
CFLAGS=$(pkg-config --cflags qmc)
LIBS=$(pkg-config --libs qmc)

echo "== pkg-config: $(pkg-config --modversion qmc)"
case "$LIBS $CFLAGS" in
  *fsanitize*) echo "FAIL: qmc.pc leaks sanitizer flags: $LIBS"; exit 1 ;;
esac

echo "== installed headers compile standalone"
fail=0
# Optional-dependency headers (HDF5/NetCDF) and *_internal.h (backend private)
# are skipped; they need third-party headers a plain consumer lacks
for h in $(cd "$PREFIX/include" && find qmc -name '*.h' | sort); do
  case "$h" in
    */hdf5_writer.h|*/netcdf_writer.h) continue ;;
  esac
  printf '#include <%s>\nint main(void){return 0;}\n' "$h" > "$TMP/hdr.c"
  if ! $CC -fsyntax-only -Wall $CFLAGS "$TMP/hdr.c" 2> "$TMP/err.txt"; then
    echo "FAIL: <$h>"; head -5 "$TMP/err.txt"; fail=1
  fi
done
[ "$fail" -eq 0 ] || exit 1

echo "== shared link"
$CC -Wall "$HERE/consumer.c" $CFLAGS $LIBS -o "$TMP/consumer_shared"
LD_LIBRARY_PATH="$PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" "$TMP/consumer_shared"

echo "== static link"
STATIC_LIBS=$(echo "$LIBS" | sed 's/-lqmc//')
$CC -Wall "$HERE/consumer.c" $CFLAGS "$PREFIX/lib/libqmc.a" $STATIC_LIBS -o "$TMP/consumer_static"
"$TMP/consumer_static"

echo "install check passed"
