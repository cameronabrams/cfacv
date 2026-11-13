# Installing

## Prerequisites

- a C compiler (gcc 11 through 16 are tested)
- `swig`
- Tcl 8.6 headers
- the GNU Scientific Library, only for the single-sweep reconstruction tool
  `Reconstruct`

## Building

```bash
git clone https://github.com/cameronabrams/cfacv.git
cd cfacv/src
make config     # shows the Tcl, GSL and swig it found
make all        # or: make cfacv.so genericdataspace.so  (no GSL needed)
```

The build puts `cfacv.so` and `genericdataspace.so` in `lib/`, and
`ForcesFromLog`, `Reconstruct` and `FEPFromForces` in `bin/`.

The makefile finds Tcl, GSL and swig itself. Override what it finds with
make variables:

| variable | meaning |
|---|---|
| `TCL_INC` | directory containing `tcl.h` |
| `GSL_PREFIX` | GSL installation prefix (with `include/` and `lib/`) |
| `SWIG` | path to the `swig` executable |
| `OPT` | optimization flags (default `-O3`) |
| `DEBUG=1` | build with `-O0 -g` |
| `PROFILE=1` | build with `-pg` |

## Checking the build

```bash
tclsh tests/cv_values.tcl
```

This loads `lib/cfacv.so` into a plain Tcl interpreter, with no NAMD
involved, and for every CV type checks the computed value against an
independent calculation and every restraining-force component against a
finite difference of the restraint energy. It prints
`cv_values: 68 passed, 0 failed` and exits 0 on success.

## NAMD and Tcl

`cfacv.so` is loaded by the Tcl interpreter inside NAMD, so build it against
the same Tcl major and minor version that NAMD carries: 8.6 for NAMD 3. A
Linux NAMD 3.0.3 multicore binary loads it as built above. If a NAMD build
cannot load it, the most reliable fix is to build NAMD against the system Tcl;
`namd/Linux-x86_64.tcl` in the repository is a NAMD arch file for that.

The string method uses NAMD's replica feature (`+partitions`), which needs
an MPI build of NAMD; a multicore build refuses to start it.
