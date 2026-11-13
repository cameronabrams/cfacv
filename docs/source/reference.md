# Reference

## NAMD configuration variables

Set these with `set` in the NAMD configuration file, before the
`tclforcesscript` line. They are read by `tcl/cfacv_tclforces.tcl`.

| variable | required | meaning |
|---|---|---|
| `CFACV_BASEDIR` | | the cfacv checkout; falls back to the `CFACV_BASEDIR` environment variable, then `$HOME/cfacv` |
| `labelPDB` | yes | the label PDB that defines the groups |
| `cvINP` | yes | the CV file |
| `restrINP` | one of these two | the restraint file |
| `restrPARAMS` | one of these two | options for one restraint per CV, used when there is no `restrINP` |
| `seed` | | random-number seed for the TAMD noise; default from the clock |
| `restartINP` | | file of starting $z$ values, one per restraint |
| `TAMDof` | | report to the log every this many steps (default 1) |
| `TAMDoutputlevel` | | what to report; a sum of the flags below (default 3) |
| `TAMDoutputFile` | | write the report to this file instead of the log |
| `TAMDbinOutputFile` | | also write a binary report to this file |
| `TAMDbinof` | | binary report every this many steps (default 1) |
| `TAMDverbose` | | if set at all, print the groups and CVs as they are set up |

(input-files)=
## Input files

cfacv reads three input files: the label PDB (`labelPDB`), the CV file
(`cvINP`) and the restraint file (`restrINP`). Their formats are a public
interface, which other programs may write. Within a major version they
change only compatibly: a file that is valid under this specification stays
valid, and means the same thing. `tests/input_files.tcl` checks cfacv
against this section.

### Label PDB

A PDB file congruent with the NAMD structure, atom for atom, that assigns
atoms to groups. A group's center is its center of mass, and every CV is
defined on group centers.

- `ATOM` and `HETATM` records are read; every other record is ignored.
- **Serial** (columns 7-11): the atom's 1-based index in the NAMD structure.
  Above 99999 it may be written in hexadecimal, as VMD writes it.
- **Occupancy** (columns 55-60): the atom's mass in amu. Only atoms in a
  group need a mass; the group's mass is the sum.
- **Beta** (columns 61-66): the atom's group, as a whole number. `0` (written
  any way, such as `0.00`, `0.0` or `-0.00`) means no group. A value `i >= 1`
  puts the atom in group `i`, which the CV file calls group `i - 1`.
- Each atom is in at most one group, so groups cannot overlap.
- The groups used must be numbered `1` to `N` with none skipped.

A beta that is not a whole number, a negative beta, a group atom without a
numeric occupancy, or a gap in the group numbers is an error.

### CV file

One CV per line; blank lines and lines starting with `#` are ignored:

```text
TYPE g1 g2 ...
```

`TYPE` is one of the {ref}`CV types <cv-types>` and `g1 g2 ...` are
**0-based** group indices: group `i` in the label PDB is `i - 1` here. The
CVs are numbered from 0 in file order.

### Restraint file

One restraint per line; blank lines and lines starting with `#` are ignored.
A restraint acts on a linear combination of the CVs: the line starts with one
coefficient per CV, in CV order, followed by the restraint's options as
`{key value}` pairs (see {ref}`restraint options <restraint-options>`):

```text
c1 c2 ... cN {key value} {key value} ...
```

For example, with two CVs, a periodic restraint on the first:

```text
1 0 {k 100.0} {rf PERIODIC} {b 1.0261928232448847}
```

Values are read as written, so targets keep their full precision. Without a
restraint file, `restrPARAMS` puts one restraint with the given options on
each CV.

(cv-types)=
## CV types

| type | groups | value | unit |
|---|---|---|---|
| `BOND` | 2 | distance between the centers | Å |
| `ANGLE` | 3 | angle at the middle center | radian |
| `DIHED` | 4 | dihedral angle, in $[-\pi, \pi]$ | radian |
| `RGYR` | any number | mass-weighted radius of gyration of the centers | Å |
| `CARTESIAN_X`, `CARTESIAN_Y`, `CARTESIAN_Z` | 1 | one coordinate of the center | Å |

The dihedral follows the IUPAC sign convention, the same as VMD's
`measure dihed`: looking along the bond from the second center to the third,
the angle is positive when the bond to the first center must rotate clockwise
to eclipse the bond to the fourth.
`RGYR` uses the group masses taken from the label PDB's occupancy column.

(restraint-options)=
## Restraint options

Each option has a long and a short key; either may be used.

| key | short | default | meaning |
|---|---|---|---|
| `SpringConstant` | `k` | 0 | spring constant $\kappa$, kcal/mol per (CV unit)² |
| `TargetValue` | `b` | 0 | target of a fixed restraint |
| `RestraintFunction` | `rf` | `HARMONIC` | `HARMONIC` or `PERIODIC` |
| `Minimum`, `Maximum` | `min`, `max` | 0 | domain bounds; set to $\pm\pi$ automatically for `PERIODIC` |
| `TAMDkT` | `kT` | | TAMD temperature $k_B\bar T$, kcal/mol |
| `TAMDgamma` | `gamma` | | TAMD friction $\bar\gamma$ (see {doc}`tamd`) |
| `TAMDdt` | `dt` | | TAMD time step, ps; normally NAMD's time step |
| `TMDTargetValue` | `tmd_z` | | final target of a targeted-MD restraint |
| `TMDInitialStep` | `tmd_t0` | | step at which the target starts to move |
| `TMDFinalStep` | `tmd_t1` | | step at which it reaches `tmd_z` |

For a restraint with value $\theta$ (the linear combination of CVs) and
target $z$, the harmonic restraint energy is $\tfrac12\kappa(\theta - z)^2$.
The periodic restraint uses the same form with $\theta - z$ wrapped into
$[-\pi, \pi]$.

## Output

Every `TAMDof` steps cfacv writes lines tagged by `TAMDoutputlevel`, a sum of
these flags:

| flag | tag | content, one column per restraint |
|---|---|---|
| 1 | `CFACV/C) Z` | the target $z$ (the TAMD variable) |
| 2 | `CFACV/C) Th` | the restraint value $\theta$ |
| 4 | `CFACV/C) FD` | the restraint force $-\kappa(\theta - z)$ |
| 8 | `CFACV/C) ND` | the noise term of the last TAMD update |

Each line gives the step number, then the values. The default, 3, writes `Z`
and `Th`. cfacv also writes `CFACV) Average gradF`, the restraint force
averaged since the last report, and `CFACV) Average MT`, the averaged metric
tensor, which the string method uses.
