#!/bin/bash
tar zvcf output.tgz init_* go*job0* restr*job0* output *fep_* image_forces.dat fep.log fep.dat fep*pdf fep.tex fep.aux rmsd.aux  rmsd-crop.pdf  rmsd.dat  rmsd-inc.pdf  rmsd.log  rmsd.pdf  rmsd.tex
rm -rf init_* go*job0* restr*job0* output *fep_* image_forces.dat fep.log fep.dat fep*pdf fep.tex fep.aux rmsd.aux  rmsd-crop.pdf  rmsd.dat  rmsd-inc.pdf  rmsd.log  rmsd.pdf rmsd.tex 



