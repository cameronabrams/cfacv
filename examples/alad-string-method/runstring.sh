#!/bin/bash
#
# String method in collective variables via NAMD/replica
#
# (c) 2016 Cameron F Abrams, Drexel University
#
CHARMRUN=/home/cfa/namd/NAMD_2.11_Source/Linux-x86_64-g++/charmrun
NAMD2=/home/cfa/namd/NAMD_2.11_Source/Linux-x86_64-g++/namd2 
NP=24  # number of processors to use for each MD simulation
NI=24
${CHARMRUN} -n $NP ${NAMD2} +replicas ${NI} job0.conf +stdout output/%d/job0.%d.log

