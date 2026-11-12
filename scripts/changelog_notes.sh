#!/usr/bin/env bash
# Print the body of one CHANGELOG.md section, without its heading.
#
# Usage: scripts/changelog_notes.sh <version|Unreleased> [CHANGELOG.md]
#   scripts/changelog_notes.sh 1.1.0     # also accepts v1.1.0
#   scripts/changelog_notes.sh Unreleased
#
# Prints nothing (and exits 1) when the section is missing.
set -euo pipefail
WANT="${1:?Usage: scripts/changelog_notes.sh <version|Unreleased> [CHANGELOG.md]}"
WANT="${WANT#v}"
FILE="${2:-CHANGELOG.md}"
awk -v want="$WANT" '
    /^## \[/ {
        if (found) exit
        v = $0; sub(/^## \[/, "", v); sub(/\].*$/, "", v)
        if (v == want) { found = 1; next }
    }
    found { lines[++n] = $0 }
    END {
        if (!found) exit 1
        # trim leading and trailing blank lines
        s = 1; while (s <= n && lines[s] ~ /^[[:space:]]*$/) s++
        e = n; while (e >= s && lines[e] ~ /^[[:space:]]*$/) e--
        for (i = s; i <= e; i++) print lines[i]
    }
' "$FILE"
