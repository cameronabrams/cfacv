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
grep -w "SM iter $iter" $log | grep -w reparam | awk '{print $9,$10}' | head -12 > s.dat
grep -w "SM iter $iter" $log |  grep -w gradient | awk '{print $9,$10}' | head -12 > gg.dat
paste -d' ' s.dat gg.dat > g.dat
../../mov/mkplt.gnu
convert -density 120 tmp.eps -flatten -crop 730x730+60+0 raw${nn}.ppm
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
 ffmpeg -r 60 -f image2 -s 730x630 -i raw%04d.ppm -vcodec libx264 -crf 25 -pix_fmt yuv420p test.mp4
endif
