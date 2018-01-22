#!/bin/tcsh
set mkmov=1
if ( $#argv > 0 ) set mkmov=1

set log=alad_sm.job0.0.history
set diter=100

@ niter = `tail -1 $log | awk '{print $3}'`
@ nf = `echo $niter / $diter | bc`

@ n = 0

cat > gtmp.gnu << EOF
TERM_W=3.35
TERM_H=2.95
set term cairolatex pdf standalone \
         color transparent crop \
         lw 2.0 rounded dl 1.5 \
         fontscale 0.7 \
         size TERM_W,TERM_H

set output "still.tex"
set xr [-pi:pi]
set yr [-pi:pi]
set zr [0:16]
set cbr [0:16]

set xlabel '\$\phi\$ (rad)'
set ylabel '\$\psi\$ (rad)' offset 0.5,0

unset key

set format cb "%.1f"

set palette rgbformulae 33,13,10
set palette maxcolors 200 

s=0.1

p "../../../alad-single-sweep/output/femap.dat" not w image,\
EOF

set iter=0
while ( $n < $nf ) 

 set nn=$n
 if ( $n < 1000 ) set nn="0"$nn
 if ( $n < 100 ) set nn="0"$nn
 if ( $n < 10 ) set nn="0"$nn

echo "SM iter $iter"
set sdat=s${iter}.dat
if ( ! -f $sdat ) then
grep -w "SM iter $iter" $log | grep -w reparam | awk '{print $9,$10}' | head -12 > s${iter}.dat
endif
set frac=`echo "scale=3; $n / $nf " |bc -l`

cat >> gtmp.gnu << EOF
    '$sdat' not w lp pt 7 ps 0.2 lc palette frac $frac, \
EOF

echo $nn

@ n++

@ iter += $diter
end
cat >> gtmp.gnu << EOF
EOF

gnuplot gtmp.gnu
pdflatex still.tex
pdfcrop still.pdf
convert -density 200 still.pdf -flatten -resize 600x500 still.png

