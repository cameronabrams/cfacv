# Changelog

All notable changes to CFACV are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[Semantic Versioning](https://semver.org/).

Releases are made as needed (see `RELEASING.md`). Add an entry under
`[Unreleased]` in the same pull request as the change it describes.

## [Unreleased]

### Fixed
- The C library builds again on current compilers: gcc 14 and later make
  implicit function declarations an error, and `strcmp` and `isspace` were
  called without their headers (#3, Prince Gyasi).
- Swapped `printf` arguments in the metric-tensor dump, and functions that
  could return no value on a guard-failure path whose callers branch on the
  result (#3).
- The `alad-string-method` example runs: `stringmethod.tcl` no longer aborts
  on an unset `nu`, hardcoded paths are gone, `mkinitimages.sh` no longer
  needs `bc`, and the restart configuration matches the protocol the driver
  implements (#3).

### Changed
- `src/makefile` finds Tcl, GSL and swig itself; `make config` reports what
  it found, and `TCL_INC`, `GSL_PREFIX` and `SWIG` override it (#3).
- The version is kept in one place, the top-level `VERSION` file, which both
  the Tcl banner and the C library report. Previously they carried separate,
  hardcoded numbers (1.00 and 0.10).
- The default branch is `main`.

### Added
- Continuous integration: every pull request builds CFACV with gcc 11 and
  gcc 14 and runs `tests/cv_values.tcl`, which checks each CV type's value
  and restraining forces without NAMD (#5).
- An MIT license, a `CITATION.cff`, this changelog, and release tooling
  (`scripts/release.sh`).
- The README is rewritten for the current tree: installation with the new
  makefile, NAMD 3 and Tcl, a table of CV types, and what each directory holds.

## [1.0.1] - 2016-03-04

Release with Zenodo integration.

## [1.0] - 2016-03-04

First public release.
