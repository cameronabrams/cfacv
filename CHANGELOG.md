# Changelog

All notable changes to cfacv are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[Semantic Versioning](https://semver.org/).

Releases are made as needed (see `RELEASING.md`). Add an entry under
`[Unreleased]` in the same pull request as the change it describes.

## [Unreleased]

### Changed
- The project name is written in lower case, cfacv, throughout the
  documentation, the citation metadata and release titles. Log tags
  (`CFACV)`, `CFACV/C)`) and `CFACV_*` variables are unchanged.
- The installation docs say which tools need GSL (only `Reconstruct`), give
  the package to install on common systems, and show how to build without it.

### Fixed
- `sbin/measurefep.sh`, which computes the free energy along a converged
  string, could not run as shipped: it hardcoded paths to a NAMD 2.11 build
  and `charmrun`, needed `bc`, ignored `-smh`, and exited with success on
  errors. It is rewritten in bash and awk: it runs `namd3 +pN` by default
  (`-namd`, `-np`, `-cr` choose otherwise), reads restarts' histories, honors
  `CFACV_BASEDIR`, resumes past images that already finished, has a dry-run
  mode (`-n`), and exits nonzero on failure. On a full 12-image run its
  restraint anchors, image forces and `fep.dat` are identical to those of the
  same steps done by hand.
- The string-method example's mean-force runs (`fepX.conf`) ran without PME:
  `set PME on` set a Tcl variable, not the NAMD option, so the mean forces
  were measured with cutoff electrostatics while the string itself was run
  with PME. `fepX.conf` now turns PME on.
- An explicit `GSL_PREFIX` now takes precedence over a `gsl-config` on the
  `PATH`; before, `gsl-config` won and `GSL_PREFIX` was ignored.
- The string method's convergence monitor, `sbin/measurermsd.sh`, which the
  example README tells you to run, did not exist. The two scripts that did,
  `measurermsd_JMG.sh` (needs `bc`) and `sbin/rmsd.sh` (could not find the
  default history file), are replaced by one `measurermsd.sh` that needs only
  bash and awk, and reads restarts' histories too.
- `FEPFromForces` no longer links GSL, which it never used, so only
  `Reconstruct` needs GSL. It now exits with a message, instead of crashing,
  when `-i` or `-o` is missing, a file cannot be opened, or the input has no
  data; its error message named the option `-d` when it is `-dim`.

## [1.1.1] - 2026-11-12

### Added
- Documentation at <https://cfacv.readthedocs.io>: installation, NAMD
  configuration variables, CV types, restraint options, output, the TAMD
  update as implemented, and the examples. CI builds it on every pull request.
- A documentation page on the string method: how it runs, dual-image mode,
  climbing, every parameter, and the limitation that a string of periodic
  CVs must not cross ±π.
- `tests/string_method.tcl`, which checks the image update, reparameterization
  against an independent implementation, fixed ends, and dual-image mode.
- `tests/tmd.tcl`, which checks the targeted-MD schedule for harmonic and
  periodic restraints, including a target reached across ±π. CI and
  `scripts/check.sh` now run every test in `tests/`.

### Fixed
- The string-method example's `runstring.sh` and `mkinitimages.sh` overwrote
  an exported `CFACV_BASEDIR`, so the example failed when staged outside a
  checkout (as on a cluster's scratch space). They now use it when set.
- String-method reparameterization assigned pointers where it meant to copy
  values, so after the first pass every reparameterization interpolated in
  place, reading images it had already moved. The resulting string was
  slightly off the reparameterized curve (by about 2e-4 after 20,000 noisy
  updates in a test of the example's settings).
- String-method reparameterization stopped the run with "did not converge"
  when it converged on exactly the last allowed iteration.
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
- Continuous integration: every pull request builds cfacv with gcc 11 and
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
