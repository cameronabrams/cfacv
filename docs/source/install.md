# Installing

## Prerequisites

- a C compiler (gcc 11 through 16 are tested)
- `swig`
- Tcl 8.6 headers
- the GNU Scientific Library (GSL) development files, only for the
  single-sweep reconstruction tool `Reconstruct`

### GSL

Only `Reconstruct` uses GSL (its QR solver and its Mersenne Twister random
number generator). `cfacv.so`, `genericdataspace.so`, `ForcesFromLog` and
`FEPFromForces` build without it, so TAMD, TMD and the string method do not
need it.

`Reconstruct` needs the headers to build and the shared library to run. A
distribution's runtime package alone (`libgsl28`, `libgsl27`, ...) is not
enough to build; install the development package:

| system | command |
|---|---|
| openSUSE | `sudo zypper install gsl-devel` |
| Fedora, RHEL, Rocky | `sudo dnf install gsl-devel` |
| Debian, Ubuntu | `sudo apt install libgsl-dev` |
| macOS (Homebrew) | `brew install gsl` |
| conda | `conda install -c conda-forge gsl` |

On a cluster, a GSL module (`module load gsl`, or similar) usually provides
`gsl-config`. The makefile looks for GSL in this order: `GSL_PREFIX` if you set
it, then `gsl-config` on the `PATH`, then the active conda environment,
then the compiler's default paths. It adds an rpath, so `Reconstruct` finds
the library at run time without `LD_LIBRARY_PATH`. `make config` reports
which GSL it found, or says it found none.

Without GSL, `make all` stops at `Reconstruct`. Build everything else by
name:

```bash
make cfacv.so genericdataspace.so ForcesFromLog FEPFromForces
```

## Building

```bash
git clone https://github.com/cameronabrams/cfacv.git
cd cfacv/src
make config     # shows the Tcl, GSL and swig it found
make all        # GSL needed for Reconstruct only; see above
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
