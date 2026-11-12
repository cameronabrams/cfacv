# The string method

CFACV implements the on-the-fly string method in collective variables
(Maragliano and Vanden-Eijnden, *Chem. Phys. Lett.* **446**, 182 (2007)). It
finds a minimum free-energy path in CV space between two states.

The path is represented by a *string* of $R$ *images* $z^0, \dots, z^{R-1}$,
points in CV space. Each image is tethered by harmonic restraints (spring
constant $\kappa$) to an atomistic MD system, a NAMD replica. Every
`steps_per_run` MD steps, CFACV gathers from each replica the mean force and
the metric tensor accumulated since the last update, moves the images, and
reparameterizes the string so the images are equally spaced again.

The worked example is `examples/alad-string-method/`, alanine dipeptide in
$(\phi, \psi)$; its README walks through setting up initial images, running,
monitoring convergence, and computing a free-energy profile along the final
string.

## Requirements

The string method runs each image as a NAMD replica (`+partitions`), so it
needs an **MPI build of NAMD**. A multicore build refuses to start with
"+partitions other than 1 is not allowed for multicore build". The number of
MPI ranks must be a multiple of the number of replicas.

## How it runs

The NAMD configuration (in the example, `job0.conf`) sets the string-method
parameters and then sources `tcl/stringmethod.tcl`. Each replica's MD
configuration (`alad_base.namd`) loads CFACV through
`tclforcesscript $CFACV_BASEDIR/tcl/cfacv_tclforces_stringmethod.tcl`, with a
restraint spring constant given by `restrPARAMS`, for example `{{k 100}}`.

Replica 0 holds the string. In each iteration:

1. every replica runs `steps_per_run` MD steps restrained at its current image,
   accumulating the restraint force and the metric tensor;
2. replica 0 collects from every replica the averaged gradient
   $\nabla F(z) \approx \langle -\kappa(\theta - z) \rangle$ and metric
   tensor $M(z)$;
3. replica 0 moves each image,
   $z \leftarrow z - h\, M(z)\, \nabla F(z)$, with $h$ = `SMPARAMS(stepsize)`;
4. replica 0 reparameterizes the string to equal arc length (below), and
   sends each replica its new image.

`steps_per_run 1` is the true on-the-fly method: the string moves at every MD
step, with $h$ small enough that the averages and the string evolve together.

### Reparameterization

The images are placed at equal arc length along the piecewise-linear curve
through the updated images, by linear interpolation. The end images are not
moved by reparameterization. The interpolation is repeated until the standard
deviation of the inter-image distances is at most `SMPARAMS(reparam_tol)`;
if that takes more than `SMPARAMS(maxiter)` passes, the run stops with an
error.

```{warning}
Arc length is measured with the ordinary Euclidean distance in CV space, and
periodic CVs are not wrapped, either in the arc length or in the image update.
A string of dihedral angles must not cross $\pm\pi$.
```

### Dual-image mode

With `SMPARAMS(dual) 1`, each image has two replicas, so `num_replicas` is
twice the number of images. Replicas $0, \dots, R-1$ are the images; replica
$i + R$ is the partner of image $i$ and is restrained at the same point. Image
$i$ is moved with its own mean force $\nabla F$ and its partner's metric
tensor $M$, so that the two estimates are statistically independent, as the
method requires.

### Climbing

With `SMPARAMS(nu)` greater than 0, the last image climbs: before each update
(from the third iteration on), the component of its previous update along the
string's tangent at that end is scaled by $-\nu$ and added to it. Use this
with `SMPARAMS(evolve_ends) 1`.

## Parameters

Set in the string-method configuration file (in the example,
`alad_stringmethod.conf`), which `job0.conf` sources.

| variable | meaning |
|---|---|
| `num_replicas` | number of NAMD replicas: the number of images, or twice that in dual-image mode |
| `steps_per_run` | MD steps between string updates; 1 is the on-the-fly method |
| `num_runs` | number of string updates |
| `runs_per_frame` | updates per output frame (trajectory and history output) |
| `frames_per_restart` | frames between restart files |
| `TEMP` | MD temperature, K |
| `SMPARAMS(stepsize)` | $h$, the image update step |
| `SMPARAMS(evolve_ends)` | 1 to move the end images, 0 to hold them fixed |
| `SMPARAMS(reparam_tol)` | reparameterization tolerance on the spread of inter-image distances |
| `SMPARAMS(maxiter)` | maximum reparameterization passes per update |
| `SMPARAMS(dual)` | 1 for dual-image mode |
| `SMPARAMS(nu)` | climbing factor for the last image; 0 (the default) disables climbing |
| `SMPARAMS(outputlevel)` | above 0, log the inter-image distances before and after each reparameterization |
| `init_restr_file` | format of the per-image restraint files that give each image's starting point, as `b` values |
| `namd_config_file` | the MD configuration common to all replicas |

## Output and restarts

Replica 0 writes a history file with, every `runs_per_frame` updates, each
image's gradient, metric tensor, and reparameterized position. Every
`frames_per_restart` frames, every replica writes NAMD restart files and
replica 0 writes a restart script and updated per-image restraint files;
`job1.conf` in the example shows how to continue from them.

To monitor convergence, run `sbin/measurermsd.sh` in the run directory. It
reads the history files and writes `rmsd.dat`: for each frame, the RMSD of the
string from the latest string,
$\sqrt{\tfrac{1}{R}\sum_i |z^i - z^i_{\mathrm{ref}}|^2}$. The string has
converged when this stops falling. `-nA N` uses the average of the last `N`
frames as the reference; `-h` lists the other options.
