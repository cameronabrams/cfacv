#!/bin/tcsh
set mkmov=1
if ( $#argv > 0 ) set mkmov=1

set log=alad_sm.job0.0.history
set diter=50

@ niter = `tail -1 $log | awk '{print $3}'`
@ nf = `echo $niter / $diter | bc`

@ n = 0

set iter=0
while ( $n < $nf ) 

 set nn=$n
 if ( $n < 1000 ) set nn="0"$nn
 if ( $n < 100 ) set nn="0"$nn
 if ( $n < 10 ) set nn="0"$nn

echo "SM iter $iter"
if ( ! -f raw${nn}.ppm ) then
grep -w "SM iter $iter" $log | grep -w reparam | awk '{print $9,$10}' | head -12 > s${iter}.dat
grep -w "SM iter $iter" $log |  grep -w gradient | awk '{print $9,$10}' | head -12 > gg${iter}.dat
paste -d' ' s${iter}.dat gg${iter}.dat > g${iter}.dat
set sdat=s${iter}.dat
gnuplot << EOF
TERM_W=3.35
TERM_H=2.95
set term cairolatex pdf standalone \
         color transparent crop \
         lw 2.0 rounded dl 1.5 \
         fontscale 0.7 \
         size TERM_W,TERM_H


set output "s.tex"
set xr [-pi:pi]
set yr [-pi:pi]
set zr [0:16]
set cbr [0:16]

set xlabel '\$\phi\$ (rad)'
set ylabel '\$\psi\$ (rad)' offset 0.5,0

unset key

set sty l 1 lc rgb "blue" pt 7 ps .15 lw 2
set sty l 2 lc rgb "yellow" pt 7 ps .2 lw 1

set format cb "%.1f"

set palette rgbformulae 33,13,10
set palette maxcolors 48 

s=0.1

p "../../../alad-single-sweep/output/femap.dat" not w image,\
  '$sdat' not w p ls 2, '$sdat' not w lp ls 1

EOF
pdflatex s.tex
pdfcrop s.pdf
convert -density 200 s.pdf -flatten -resize 568x500 raw${nn}.ppm

endif

#convert untitled.0${nn}.ppm -resize 60x60 inset0.ppm
#foreach d (1 2 3 4 5 6 7)
#   convert ../${d}/untitled.0${nn}.ppm -resize 60x60 ../${d}/inset0.ppm
#end

#convert raw${nn}.ppm -draw 'image over 460,40 60,60 "inset0.ppm"' \
#                     -draw 'image over 460,100 60,60 "../1/inset0.ppm"' \
#                     -draw 'image over 460,160 60,60 "../2/inset0.ppm"' \
#                     -draw 'image over 460,220 60,60 "../3/inset0.ppm"' \
#                     -draw 'image over 460,280 60,60 "../4/inset0.ppm"' \
#                     -draw 'image over 460,340 60,60 "../5/inset0.ppm"' \
#                     -draw 'image over 460,400 60,60 "../6/inset0.ppm"' \
#                     -draw 'image over 460,460 60,60 "../7/inset0.ppm"' \
#                     ${nn}.ppm

echo $nn

@ n++

@ iter += $diter
end
#

if ( $mkmov ) then 
 if ( -f test.mpg ) then
   rm test.mpg
 endif
 ffmpeg -r 60 -f image2 -s 568x500 -i raw%04d.ppm -vcodec libx264 -crf 25 -pix_fmt yuv420p test.mp4
endif
