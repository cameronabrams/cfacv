# Releasing cfacv

cfacv is archived at Zenodo, which mints a DOI for every GitHub Release of
this repository. A release is therefore a published, citable artifact, and is
cut deliberately.

## When

As needed: when `[Unreleased]` in `CHANGELOG.md` holds changes worth
publishing. There is no fixed schedule.

## Day to day

Every change arrives as a pull request, and `main` accepts it only when CI
passes (`.github/workflows/tests.yml`: a build with gcc 11 and gcc 14, and
`tests/cv_values.tcl`). Add a `CHANGELOG.md` entry under `[Unreleased]` in
the same pull request.

## Cutting a release

From a clean `main` that is level with `origin/main`:

    scripts/release.sh 1.1.0

The script shows the pending entries and asks you to confirm the number,
builds the committed tree and runs the tests (`scripts/check.sh`), then
rotates `CHANGELOG.md`, writes the version to `VERSION` and `CITATION.cff`,
commits, tags `v1.1.0`, and pushes. Pushing the release commit to the
protected `main` needs a repository admin.

The pushed tag runs `.github/workflows/release.yml`, which rebuilds and
retests the tag, checks that it matches `VERSION`, and creates the GitHub
Release from the changelog section. Zenodo picks the release up from there.

## Versions

Semantic versioning, `v`-prefixed tags (`v1.1.0`). The two releases from
2016 are tagged `1.0` and `1.0.1` and are left as they are.

- MAJOR: a NAMD-config key, the `cv.inp` or `restr.inp` syntax, or a Tcl
  procedure that a config calls changed incompatibly.
- MINOR: new capability, such as a new CV type.
- PATCH: fixes only.
