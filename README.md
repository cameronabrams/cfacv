# CFACV -- a collective variables implementation for NAMD

[![tests](https://github.com/cameronabrams/cfacv/actions/workflows/tests.yml/badge.svg)](https://github.com/cameronabrams/cfacv/actions/workflows/tests.yml)
[![DOI](https://zenodo.org/badge/20736/cameronabrams/cfacv.svg)](https://zenodo.org/badge/latestdoi/20736/cameronabrams/cfacv)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Cameron F Abrams, Drexel University, Philadelphia, Pennsylvania  
cfa22@drexel.edu  
215 895 2231  

CFACV is a simple collective variables implementation for NAMD that uses
`tclforces` and implements Temperature-Accelerated Molecular Dynamics (TAMD).
Collective variables (CVs) are functions of atomic Cartesian coordinates;
CFACV's CVs are built on centers of mass of groups of atoms, and the lengths,
angles and dihedrals constructed from them. CVs can be harmonically
restrained or driven with TAMD. CFACV also implements the on-the-fly string
method in CVs and single-sweep free-energy reconstruction.

## Citing CFACV

CFACV is archived at Zenodo; cite the version you used (the DOI badge above
leads to all versions, and GitHub's "Cite this repository" reads
[`CITATION.cff`](CITATION.cff)). Please also cite

> Cameron F. Abrams and Eric Vanden-Eijnden, "Large-scale conformational
> sampling of proteins using temperature-accelerated molecular dynamics,"
> *Proc. Natl. Acad. Sci. USA* **107**, 4961-4966 (2010).
> [doi:10.1073/pnas.0914540107](https://doi.org/10.1073/pnas.0914540107)

## Installing

Prerequisites: a C compiler, `swig`, Tcl 8.6 headers, and (for the
single-sweep tools only) the GNU Scientific Library.

```bash
git clone https://github.com/cameronabrams/cfacv.git
cd cfacv/src
make config     # shows the Tcl, GSL and swig it found
make all        # or: make cfacv.so genericdataspace.so  (no GSL needed)
```

The makefile finds Tcl, GSL and swig itself. If it picks the wrong one, or
none, override it: `make all TCL_INC=/dir/with/tcl.h GSL_PREFIX=/gsl/prefix
SWIG=/path/to/swig`. `DEBUG=1` builds without optimization and with symbols.

The build puts `cfacv.so` and `libgenericdataspace.so` in `lib/` and the
single-sweep tools (`ForcesFromLog`, `Reconstruct`, `FEPFromForces`) in
`bin/`.

To check the build without NAMD:

```bash
tclsh tests/cv_values.tcl
```

This loads `lib/cfacv.so` into a plain Tcl interpreter and checks the value
and restraining forces of every CV type against independent calculations.

### NAMD and Tcl

`cfacv.so` is loaded by the Tcl interpreter inside NAMD, so it should be
built against the same Tcl major/minor version that NAMD carries (8.6 for
NAMD 3). A Linux NAMD 3.0.3 multicore binary loads it as built above
(checked 2026-11), and the string-method example has run under an MPI build
of NAMD 3. If a
NAMD build cannot load it, the most reliable fix is to build NAMD against the
system Tcl; [`namd/Linux-x86_64.tcl`](namd/Linux-x86_64.tcl) is an arch file
for doing that.

The string-method example uses NAMD's replica feature (`+partitions`), which
requires an MPI build of NAMD; a multicore build refuses to start it.

## Using CFACV in a NAMD run

CFACV hooks into NAMD through `tclforces`. In the NAMD configuration file:

```tcl
set CFACV_BASEDIR    /path/to/cfacv     ;# or set the CFACV_BASEDIR environment variable
tclforces            on
tclforcesscript      $CFACV_BASEDIR/tcl/cfacv_tclforces.tcl
set labelPDB         label.pdb
set cvINP            cv.inp
set restrINP         restr.inp          ;# or set restrPARAMS; see below
set TAMDof           100                ;# report every 100 steps
```

`cfacv_tclforces.tcl` sets up the CVs and defines the `calcforces` procedure
NAMD calls at every step; the procedures it calls are in `tcl/cfacv.tcl`, and
the computation is in `src/cfacv.c`, reached through the SWIG bindings in
`lib/cfacv.so`. To understand the implementation, trace the code starting
from `cfacv_tclforces.tcl`.

CFACV needs three inputs:

* **`labelPDB`** -- a PDB congruent with the simulation system in which the
  beta field of each atom gives its group ("center"); 0 means the atom is in
  no group. You create this file; `examples/hiv1-ma/mk_tPDB.tcl` shows how to
  cluster residues into subdomain groups automatically. A group can be a
  single atom.

  **Note:** beta-field group IDs start at 1, but `cv.inp` numbers groups from
  0. Atoms with beta 1 form group 0 in `cv.inp`, beta 2 is group 1, and so on.

* **`cvINP`** -- one CV per line: a type and the groups it uses.

  | type | groups | value |
  |---|---|---|
  | `BOND` | 2 | distance |
  | `ANGLE` | 3 | angle at the middle group, radians |
  | `DIHED` | 4 | dihedral, radians, in [-π, π] (IUPAC sign, as VMD's `measure dihed`) |
  | `RGYR` | any number | mass-weighted radius of gyration of the groups |
  | `CARTESIAN_X`, `_Y`, `_Z` | 1 | one coordinate of the group's center |

  For example, the backbone φ and ψ of alanine dipeptide
  (`examples/alad-tamd/cv.inp`):

  ```
  DIHED 0 1 2 3
  DIHED 1 2 3 4
  ```

* **`restrINP`** (optional) -- one restraint per line: the coefficients of a
  linear combination of the CVs, then its options:

  ```
  1 0 {k 100.0} {rf PERIODIC} {TAMDkT 9.5} {TAMDgamma 5.0} {TAMDdt 0.001}
  0 1 {k 100.0} {rf PERIODIC} {TAMDkT 9.5} {TAMDgamma 5.0} {TAMDdt 0.001}
  ```

  Without `restrINP`, every CV gets its own restraint with the same options,
  given in the NAMD configuration file:

  ```tcl
  set restrPARAMS {{k 100} {TAMDkT 6.0} {TAMDgamma 100.0} {TAMDdt 0.002}}
  ```

  Units: `k` in kcal/mol per CV unit squared, `TAMDkT` in kcal/mol,
  `TAMDgamma` in 1/ps, `TAMDdt` in ps. A restraint on a dihedral must be
  declared periodic with `{rf PERIODIC}`.

The log reports CVs and TAMD variables on lines that begin with `CFACV)`.

## Examples

Each directory in [`examples/`](examples) is a complete input set:

* `alad-fixed-phi-psi/` -- alanine dipeptide (CHARMM22, TIP3P) restrained at
  a fixed (φ, ψ).
* `alad-tamd/` -- TAMD on alanine dipeptide's φ and ψ.
* `alad-single-sweep/` -- single-sweep reconstruction of the (φ, ψ) free
  energy surface; see its README.
* `alad-string-method/` -- the on-the-fly string method in CVs; see its
  README. Needs an MPI build of NAMD.
* `hiv1-ma/` -- TAMD on the centers of five subdomains of HIV-1 matrix
  protein (PDB 1HIW). `mk_tPDB.tcl` builds the `label.pdb` (run it with
  `vmd -dispdev text -e mk_tPDB.tcl -args -p my_system.pdb`).

To run one, `cd` into its directory and run NAMD on its configuration file,
for example `namd3 +p4 tamd.conf`. The configurations set
`CFACV_BASEDIR ../..`, so they run in place from a checkout.

## What else is here

* `tcl/` -- the tclforces drivers: `cfacv_tclforces.tcl` (TAMD and
  restraints), `cfacv_tclforces_stringmethod.tcl` with `stringmethod.tcl`,
  and `cfacv_tclforces_mil2.tcl` (milestoning; needs the PVRW patch below).
* `namd-patch/` -- the position-velocity-rewind (PVRW) patch to NAMD 2.11
  that the milestoning driver uses, with a test case. It has not been ported
  to NAMD 3.
* `namd/` -- a NAMD arch file for building NAMD against the system Tcl.
* `sbin/` -- shell helpers for post-processing runs.

## Development and releases

Changes arrive as pull requests; CI builds with gcc 11 and gcc 14 and runs
`tests/cv_values.tcl`. Record each change in [`CHANGELOG.md`](CHANGELOG.md)
under `[Unreleased]`. [`RELEASING.md`](RELEASING.md) describes how a release
is cut.

Bug reports and questions are welcome as GitHub issues.

## Acknowledgments

1. VMD and NAMD are products of the [Theoretical and Computational Biophysics
   Group at the NIH Center for Macromolecular Modeling and Bioinformatics at
   the University of Illinois at Urbana-Champaign](http://www.ks.uiuc.edu).

2. All code and data in this repository have been made possible with partial
   support from NIH through grants AI084117, AI093248, GM115249, GM056550,
   and GM100472, the National Science Foundation through grants DMR-1207389
   and MCB-1330205, and the US Army through grants W911NF-12-2-0022,
   W911-NF-13-1-0046, and W911NF-12-R-0011.

Contributors: Cameron F Abrams, Jasmine Gardner, Gourav Shrivastav, Prince
Gyasi.
