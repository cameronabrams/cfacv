#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

rmsd=rmsd.dat
   
# PLOT PES (IT USE PDFLATEX)      
gnuplot << HEREGNUPLOT

TERM_W=3.35 #8.5
TERM_H=2.95 #7.5

set term cairolatex pdf standalone \
         color \
         transparent \
         crop \
         lw 3.5 rounded dl 1.5 \
         fontscale 0.7 \
         size TERM_W,TERM_H

C0='#000000' # black
C1='#00CAFF' # blue
C2='#ED3221' # red    
C3='#FECE23' # yellow 
C4='#35E90E' # green 

set output "rmsd.tex"
set ylabel 'CV-distance'
set xlabel 'iteration number'
set format y '%.1f'
set xr [0:*]
set yr [0:*]

plot '$rmsd' u 1:2  w l lc rgb C1 notit

HEREGNUPLOT

pdflatex rmsd.tex
pdfcrop rmsd.pdf
