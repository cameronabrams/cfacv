# Build notes

Written 2026-11-05. Not upstream; this file is local to the `pkgyasi` fork and
is untracked unless you choose to commit it.

## Two builds exist

| where | path | compiler | status |
|---|---|---|---|
| picotte | `/ifs/groups/abramsGrp/pg565/cfacv` | gcc 9.2.0 | builds clean, **this is the production one** |
| this workstation | `~/softwares/cfacv` | gcc 15.3.0 | **builds clean since 2026-11-09**, see below |

Both produce `lib/cfacv.so`, `lib/genericdataspace.so` and the three binaries
`bin/ForcesFromLog`, `bin/Reconstruct`, `bin/FEPFromForces`. Both were verified by loading the shared objects into `tclsh`.

## Loading it is not what the filename suggests

    load lib/cfacv.so cfa_cvlibc
    load lib/genericdataspace.so genericdataspace

The SWIG module in `src/cfacv.i` is named **`cfa_cvlibc`**, not `cfacv`, so the
exported init symbol is `Cfa_cvlibc_Init`. A bare `load lib/cfacv.so` makes Tcl
look for `Cfacv_Init` and fail with "undefined symbol". The package name is not
optional. `tcl/cfacv.tcl` line 10 already does this correctly; the trap is only
for anyone testing the object by hand.

## picotte

    module load gsl/gcc/2.7.1
    cd /ifs/groups/abramsGrp/pg565/cfacv && mkdir -p lib bin
    cd src && make all

Nothing else is needed. picotte carries `swig` 3.0.12 system-wide and
`tcl-devel` 8.6.8 with `tcl.h` in `/usr/include`, which is where the makefile's
hardcoded `-I/usr/include/tcl` effectively resolves.

## This workstation

**As of 2026-11-09 the long flag incantation below is no longer needed.**
`src/makefile` now discovers Tcl, GSL and swig itself, and commit 4f8674c added
the missing `<string.h>` and `<ctype.h>` includes, so gcc 15.3.0 builds the tree
clean with no warning suppression:

    conda create -y -n cfacv-build swig gsl tk tcl -c conda-forge   # once
    E=~/miniconda3/envs/cfacv-build
    cd ~/softwares/cfacv/src
    PATH="$E/bin:$PATH" CONDA_PREFIX="$E" make all

`make config` prints what was discovered before anything is compiled. Verified
2026-11-09 on gcc 15.3.0: no errors, no implicit-declaration warnings, all three
binaries resolve their shared libraries, and `lib/cfacv.so` registers 171
commands including 31 `DataSpace_*` with `NewDataSpace` constructing an object.

Useful overrides:

| want | pass |
|---|---|
| tune to this cpu (not portable) | `OPT='-O3 -march=native'` |
| profile with gprof | `PROFILE=1` |
| debug build | `DEBUG=1` |
| explicit paths | `TCL_INC=/dir SWIG=/path GSL_PREFIX=/prefix` |

**`-ffast-math` is deliberately absent and should stay absent.** This code
returns forces to NAMD; relaxing IEEE semantics can change a trajectory, so it
must not be enabled without revalidating against a reference run.

The pre-2026-11 makefile is kept as `src/makefile.orig`, untracked.

### What the old incantation was, and why it is recorded

Before the makefile rework the build needed every path supplied by hand plus
`-Wno-implicit-function-declaration -Wno-int-conversion`, because the sources
called `strcmp` and `isspace` without including their headers. Under gcc 9 that
is a warning; **from gcc 14 an implicit function declaration is an error**, so
the build stopped at `cfacv.o`. The includes fix is the real remedy and is in
this fork; the flags are no longer required anywhere.

## NAMD version

`namd-patch/cfa_pvrw_namd211.patch` targets NAMD 2.11. This project runs NAMD 3,
and `/ifs/groups/abramsGrp/pg565/_cfacv_namd3_test` is a NAMD 3.0.3 probe case.
The plain TAMD path goes through `tclforces` and should not need the patch. If
the string method does need it, applying it to a NAMD 3 tree is real work and
should be scoped before Aim ii depends on it.
