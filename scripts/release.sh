#!/usr/bin/env bash
# Release CFACV at a given version.
#
# Usage: scripts/release.sh <version> [--yes]
# Example: scripts/release.sh 1.1.0
#
# Prerequisites (checked automatically):
#   - Working tree clean, on main, level with origin/main
#   - CHANGELOG.md has a non-empty "## [Unreleased]" section
#   - The tag v<version> does not already exist, locally or on origin
#   - The committed tree builds and tests/cv_values.tcl passes (scripts/check.sh)
#
# What it does:
#   1. Rotates CHANGELOG.md: [Unreleased] -> [<version>] - <date>, with a fresh
#      empty [Unreleased] above it
#   2. Writes <version> to VERSION and to version/date-released in CITATION.cff
#   3. Commits those as "Release v<version>" and tags it v<version>
#   4. Pushes the commit and the tag to origin
#
# The pushed tag triggers .github/workflows/release.yml, which rebuilds and
# retests the tagged commit and creates the GitHub Release. Zenodo archives
# every GitHub Release and mints a DOI for it, so PUSHING THE TAG PUBLISHES A
# CITABLE ARTIFACT. That is a human's decision, not an agent's.
#
# main is protected (changes arrive by pull request with CI passing). The
# release commit is pushed directly, which only a repository admin can do;
# GitHub will report that a rule was bypassed. That is expected.

set -euo pipefail

ASSUME_YES=0
ARGS=()
for arg in "$@"; do
    case "$arg" in
        --yes|-y) ASSUME_YES=1 ;;
        *) ARGS+=("$arg") ;;
    esac
done
set -- ${ARGS[@]+"${ARGS[@]}"}

VERSION="${1:?Usage: scripts/release.sh <version> [--yes]  (e.g. 1.1.0)}"
VERSION="${VERSION#v}"
TODAY="$(date +%Y-%m-%d)"

if ! [[ "$VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+([-.][0-9A-Za-z.]+)?$ ]]; then
    echo "ERROR: '$VERSION' is not a semantic version (e.g. 1.1.0)" >&2
    exit 1
fi

cd "$(git rev-parse --show-toplevel)"

# ── Preconditions ─────────────────────────────────────────────────────────────

if ! git diff --quiet || ! git diff --cached --quiet; then
    echo "ERROR: working tree has uncommitted changes; commit or stash them first" >&2
    exit 1
fi

BRANCH="$(git branch --show-current)"
if [ "$BRANCH" != "main" ]; then
    echo "ERROR: must be on main (currently on '$BRANCH')" >&2
    exit 1
fi

git fetch --quiet origin main
if [ "$(git rev-parse HEAD)" != "$(git rev-parse origin/main)" ]; then
    echo "ERROR: main and origin/main differ; pull first, so the tag names a" >&2
    echo "       commit that exists on origin" >&2
    exit 1
fi

if git rev-parse -q --verify "refs/tags/v$VERSION" >/dev/null; then
    echo "ERROR: tag v$VERSION already exists locally" >&2
    exit 1
fi
if git ls-remote --tags origin "refs/tags/v$VERSION" | grep -q .; then
    echo "ERROR: tag v$VERSION already exists on origin" >&2
    exit 1
fi

NOTES="$(scripts/changelog_notes.sh Unreleased || true)"
if [ -z "$(printf '%s' "$NOTES" | tr -d '[:space:]')" ]; then
    echo "ERROR: '## [Unreleased]' in CHANGELOG.md is missing or empty; write the" >&2
    echo "       notes for this release first. (An empty section makes a release" >&2
    echo "       whose notes read as 'nothing changed'.)" >&2
    exit 1
fi

# ── Confirm the number against what is pending ────────────────────────────────

LAST="$(git tag --list 'v*' --sort=-v:refname | head -1)"
LAST="${LAST:-$(cat VERSION)}"
echo
echo "  Pending for v$VERSION (last release: ${LAST#v}):"
echo "  ------------------------------------------------------------"
printf '%s\n' "$NOTES" | sed 's/^/  /'
echo "  ------------------------------------------------------------"
echo "  MAJOR if a config key, cv.inp/restr.inp syntax, or a Tcl procedure"
echo "  that a NAMD config calls changed incompatibly. MINOR for new"
echo "  capability. PATCH for fixes only."
echo
if [ "$ASSUME_YES" != "1" ]; then
    if [ ! -t 0 ]; then
        echo "ERROR: not a terminal, so the number cannot be confirmed; re-run with" >&2
        echo "       --yes once you have checked the entries above." >&2
        exit 1
    fi
    printf "  Is %s the right number for these changes? [y/N] " "$VERSION"
    read -r reply
    case "$reply" in
        y|Y|yes|YES) ;;
        *) echo "  Stopped; nothing has been changed." >&2; exit 1 ;;
    esac
fi

# ── Build and test what is committed ──────────────────────────────────────────
# CI does this again on the tag, but only AFTER the tag is pushed, and a tag
# that Zenodo has archived cannot be taken back.
#   - Skip (you have decided a failure is not a blocker): SKIP_TESTS=1
if [ "${SKIP_TESTS:-0}" != "1" ]; then
    echo "Building and testing (skip with SKIP_TESTS=1)..."
    if ! scripts/check.sh HEAD; then
        echo "ERROR: build or tests failed; not releasing." >&2
        exit 1
    fi
fi

# ── Edit CHANGELOG.md, VERSION, CITATION.cff ──────────────────────────────────

restore() { git checkout -- CHANGELOG.md VERSION CITATION.cff; }

if [ "$(grep -c '^## \[Unreleased\]' CHANGELOG.md)" != "1" ]; then
    echo "ERROR: expected exactly one '## [Unreleased]' in CHANGELOG.md" >&2
    exit 1
fi
sed -i "s/^## \[Unreleased\]\$/## [Unreleased]\n\n## [$VERSION] - $TODAY/" CHANGELOG.md
printf '%s\n' "$VERSION" > VERSION
sed -i -e "s/^version: .*/version: $VERSION/" \
       -e "s/^date-released: .*/date-released: $TODAY/" CITATION.cff

# Every edit must have landed, exactly once, before anything is committed.
if [ "$(scripts/changelog_notes.sh "$VERSION")" != "$NOTES" ] \
   || [ "$(cat VERSION)" != "$VERSION" ] \
   || [ "$(grep -c "^version: $VERSION\$" CITATION.cff)" != "1" ] \
   || [ "$(grep -c "^date-released: $TODAY\$" CITATION.cff)" != "1" ]; then
    echo "ERROR: CHANGELOG.md, VERSION and CITATION.cff do not agree after the" >&2
    echo "       edit; restored them. Nothing has been committed." >&2
    restore
    exit 1
fi

# ── Commit, tag, push ─────────────────────────────────────────────────────────

git add CHANGELOG.md VERSION CITATION.cff
git commit -q -m "Release v$VERSION"
git tag -a "v$VERSION" -m "CFACV v$VERSION"

echo "Pushing the release commit and tag v$VERSION to origin..."
git push origin main
git push origin "v$VERSION"

echo
echo "Done. release.yml will rebuild, retest, and create the GitHub Release;"
echo "Zenodo then archives it and mints the DOI for this version."
echo "  https://github.com/cameronabrams/cfacv/releases/tag/v$VERSION"
