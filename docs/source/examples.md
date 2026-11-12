# Examples

Each directory under `examples/` in the repository is a complete input set.
To run one, `cd` into it and run NAMD on its configuration file, for example
`namd3 +p4 tamd.conf`. The configurations set `CFACV_BASEDIR` to `../..`, so
they run in place from a checkout.

| directory | what it shows | configuration |
|---|---|---|
| `alad-fixed-phi-psi/` | alanine dipeptide (CHARMM22, TIP3P) restrained at a fixed $(\phi, \psi)$ | `fixed-phi-psi.conf` |
| `alad-tamd/` | TAMD on alanine dipeptide's $\phi$ and $\psi$ | `tamd.conf` |
| `alad-single-sweep/` | single-sweep reconstruction of the $(\phi, \psi)$ free-energy surface; has its own README | `tamd.conf`, `fmX.conf` |
| `alad-string-method/` | the on-the-fly string method in CVs (see {doc}`string_method`); has its own README; needs an MPI build of NAMD | `alad_stringmethod.conf` |
| `hiv1-ma/` | TAMD on the centers of five subdomains of HIV-1 matrix protein (PDB 1HIW) | `example.conf` |

For `hiv1-ma/`, first build the label PDB with VMD:

```bash
vmd -dispdev text -e mk_tPDB.tcl -args -p my_system.pdb
```
