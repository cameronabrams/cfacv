# Using CFACV in a NAMD run

CFACV hooks into NAMD through `tclforces`. A minimal stanza in the NAMD
configuration file:

```tcl
set CFACV_BASEDIR    /path/to/cfacv     ;# or the CFACV_BASEDIR environment variable
tclforces            on
tclforcesscript      $CFACV_BASEDIR/tcl/cfacv_tclforces.tcl
set labelPDB         label.pdb
set cvINP            cv.inp
set restrINP         restr.inp          ;# or set restrPARAMS
set TAMDof           100
```

`cfacv_tclforces.tcl` reads the inputs, sets up the centers, CVs and
restraints, and defines the `calcforces` procedure that NAMD calls at every
step. At each step `calcforces`

1. gets the center-of-mass position of each group from NAMD,
2. computes every CV and its gradient with respect to the centers,
3. computes each restraint's value, force and energy, and advances the TAMD
   or TMD variables,
4. hands the force on each center back to NAMD as a force on that group
   (`addforce`).

Every variable the script reads is listed in {doc}`reference`.

## The label PDB (`labelPDB`)

A PDB file congruent with the simulation system (same atoms, same order)
that assigns atoms to groups:

- **beta** (columns 61-66) is the atom's group number, starting at **1**;
  `0.00` means the atom is in no group. Number groups consecutively, 1, 2,
  3, ... with no gaps.
- **occupancy** (columns 55-60) must hold the atom's **mass**. CFACV sums it
  into the group's mass.
- Only `ATOM` records are read. Atoms written as `HETATM` (often ions,
  water and ligands) cannot be grouped.

A group can be a single atom. `examples/hiv1-ma/mk_tPDB.tcl` builds a label
PDB with VMD by clustering residues into subdomains.

```{important}
Group numbers in the label PDB start at 1, but `cv.inp` numbers groups
from **0**. Atoms with beta 1 form group 0 in `cv.inp`, beta 2 is group 1,
and so on.
```

## The CV file (`cvINP`)

One CV per line: its type, then the groups it uses. Lines beginning with `#`
are comments.

```
# alanine dipeptide backbone phi and psi
DIHED 0 1 2 3
DIHED 1 2 3 4
```

The CV types are listed in {ref}`cv-types`.

## Restraints (`restrINP` or `restrPARAMS`)

A restraint acts on a linear combination of the CVs. In a restraint file,
each line gives one coefficient per CV, then the restraint's options as
`{key value}` pairs:

```
1 0 {k 100.0} {rf PERIODIC} {TAMDkT 9.5} {TAMDgamma 5.0} {TAMDdt 0.001}
0 1 {k 100.0} {rf PERIODIC} {TAMDkT 9.5} {TAMDgamma 5.0} {TAMDdt 0.001}
```

This file puts one restraint on each of two CVs. A row such as
`1 -1 {k 50}` would restrain their difference.

Without a restraint file, set `restrPARAMS` in the NAMD configuration
instead; every CV then gets its own restraint with those options:

```tcl
set restrPARAMS {{k 100} {TAMDkT 6.0} {TAMDgamma 100.0} {TAMDdt 0.002}}
```

A restraint with all three of `TAMDkT`, `TAMDgamma` and `TAMDdt` is driven by
TAMD (see {doc}`tamd`). One with all three TMD options is targeted MD: its
target starts at the restraint's value on the first step and moves by equal
increments from step `tmd_t0` to step `tmd_t1`, ending at `tmd_z`. Otherwise it is a fixed restraint at its
target value `b`. All options are in {ref}`restraint-options`.

A restraint on a dihedral must be periodic, `{rf PERIODIC}`, so that the
restraint takes the short way around the circle.

## Restarting

To continue a TAMD run from where the auxiliary variables $z$ left off, put
their values, one number per restraint and separated by spaces, in a file and
name it with `set restartINP <file>`. Without it, each $z$ starts at the
value of its restraint's CV combination at the first step.
