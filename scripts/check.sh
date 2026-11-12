#!/usr/bin/env bash
# Build the committed tree in a scratch directory and run every tests/*.tcl.
#
# Usage: scripts/check.sh [<git-ref>]      (default: HEAD)
#
# Builds from `git archive`, so it tests what is committed, never stray files
# in the working tree, and leaves no build products in the repository. Builds
# only the two Tcl libraries: the GSL tools are not needed by the test, and CI
# builds those on every pull request.
set -euo pipefail
REF="${1:-HEAD}"
cd "$(git rev-parse --show-toplevel)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/cfacv-check.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
git archive "$REF" | tar -x -C "$WORK"
make -C "$WORK/src" cfacv.so genericdataspace.so > "$WORK/build.log" 2>&1 || {
    tail -30 "$WORK/build.log" >&2
    echo "ERROR: build failed" >&2
    exit 1
}
status=0
for t in "$WORK"/tests/*.tcl; do
    out="$(tclsh "$t" "$WORK" 2>&1)" || status=1
    printf '%s\n' "$out" | grep -E '^(FAIL|[a-z_]+: [0-9]+ passed)' || true
done
exit $status
