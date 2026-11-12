# CFACV

CFACV is a collective variables implementation for
[NAMD](https://www.ks.uiuc.edu/Research/namd/) that works through NAMD's
`tclforces` interface. It implements temperature-accelerated molecular
dynamics (TAMD), harmonic and periodic restraints on collective variables
(CVs) and linear combinations of them, targeted MD, the on-the-fly string
method in CVs, and single-sweep free-energy reconstruction.

CVs are built on *centers*, the centers of mass of groups of atoms. A CV is
a distance, angle or dihedral between centers, a radius of gyration of a set
of centers, or one Cartesian coordinate of a center.

```{toctree}
:maxdepth: 2

install
usage
reference
tamd
examples
development
```

## Citing CFACV

CFACV is archived at Zenodo; cite the version you used (GitHub's "Cite this
repository" reads the repository's `CITATION.cff`). Please also cite

> Cameron F. Abrams and Eric Vanden-Eijnden, "Large-scale conformational
> sampling of proteins using temperature-accelerated molecular dynamics,"
> *Proc. Natl. Acad. Sci. USA* **107**, 4961-4966 (2010).
> [doi:10.1073/pnas.0914540107](https://doi.org/10.1073/pnas.0914540107)

Source: <https://github.com/cameronabrams/cfacv>. CFACV is released under
the MIT license.
