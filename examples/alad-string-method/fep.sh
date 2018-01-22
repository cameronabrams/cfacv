#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

pmf=fep.dat
   
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

set output "fep.tex"
set ylabel 'free energy (kcal/mol)'
set xlabel 'CV-distance'
set format y '%.1f'
set format x '%.1f'
set xr [0:5]
set yr [0:10]

plot '$pmf' u 2:3  w l lc rgb C1 notit

HEREGNUPLOT

pdflatex fep.tex
pdfcrop fep.pdf
