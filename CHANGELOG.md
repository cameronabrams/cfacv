# Changelog

All notable changes to CFACV are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[Semantic Versioning](https://semver.org/).

Releases are made as needed (see `RELEASING.md`). Add an entry under
`[Unreleased]` in the same pull request as the change it describes.

## [Unreleased]

### Added
- Documentation at <https://cfacv.readthedocs.io>: installation, NAMD
  configuration variables, CV types, restraint options, output, the TAMD
  update as implemented, and the examples. CI builds it on every pull request.
- `tests/tmd.tcl`, which checks the targeted-MD schedule for harmonic and
  periodic restraints, including a target reached across ±π. CI and
  `scripts/check.sh` now run every test in `tests/`.

### Fixed
- Targeted MD on a periodic restraint (a dihedral) moved its target at every
  step, before `tmd_t0` and after `tmd_t1` as well as between them, and left
  it outside [-π, π] when it moved outside the window. The target now moves
  only inside the window, as it already did for harmonic restraints.
- Targeted MD applied one increment too many, so the target finished
  `(tmd_z - start)/(tmd_t1 - tmd_t0)` past `tmd_z`. It now arrives at `tmd_z`
  exactly at step `tmd_t1`. This affects harmonic and periodic restraints.
- The units of `TAMDgamma` were documented as 1/ps. The TAMD update divides
  a force by it, so it is a friction, in kcal·ps/mol per CV unit squared.
- `make config` no longer reports GSL as if found when it is not; it says
  whether GSL was found, and where.
- The README named the second library `libgenericdataspace.so`; the build
  makes, and `tcl/genericdataspace.tcl` loads, `genericdataspace.so`.

## [1.1.0] - 2026-11-12

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
