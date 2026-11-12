# TAMD as implemented

In temperature-accelerated MD, each restrained CV combination $\theta(x)$ is
coupled by a spring to an auxiliary variable $z$, and $z$ evolves by
overdamped Langevin dynamics at a temperature $\bar T$ higher than the
physical one. The atoms feel the restraint force; $z$ feels the opposite
force, and explores the free-energy surface of $\theta$ at the accelerated
temperature (Abrams and Vanden-Eijnden, *PNAS* **107**, 4961 (2010)).

At every MD step, for each TAMD restraint, CFACV

1. computes the restraint force on $\theta$,
   $f = -\kappa\,(\theta - z)$, and applies $f\,\partial\theta/\partial
   \mathbf{R}_i$ to each center $i$;
2. advances $z$ by one Euler-Maruyama step,

   $$
   z \leftarrow z + \frac{\Delta t}{\bar\gamma}\,\kappa\,(\theta - z)
              + \sqrt{\frac{2\,k_B\bar T\,\Delta t}{\bar\gamma}}\;\xi ,
   $$

   where $\Delta t$ is `TAMDdt`, $\bar\gamma$ is `TAMDgamma`, $k_B\bar T$ is
   `TAMDkT`, and $\xi$ is a unit-variance random number (the sum of six
   uniform deviates, shifted and scaled). For a periodic restraint the new
   $z$ is wrapped into $[-\pi, \pi]$.

On the first step, $z$ is set to the current $\theta$ unless a `restartINP`
gives it.

## Units of the friction

Because the update divides a force by $\bar\gamma$, $\bar\gamma$ carries the
units of a friction coefficient, not a rate: kcal·ps/mol per (CV unit)², so
that $\Delta t\,f/\bar\gamma$ is in CV units. Documentation before v1.2
labeled `TAMDgamma` as 1/ps; that reading holds only if the auxiliary
variable's mass is taken as 1 in those units. What matters for a run is the
ratio $\bar\gamma/\kappa$, the relaxation time of $z$ toward $\theta$, which
should be long compared with the time the atoms need to follow $z$.
