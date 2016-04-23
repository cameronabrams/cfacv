#!/usr/bin/gnuplot
set contour
set table "ref_md.contours"
set cntrparam levels discrete 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 6, 8, 10
unset surf
splot "ref_md.dist" u 1:2:7
unset table
unset contour

set term post eps enh color "Helvetica" 28 solid
set out "femap.eps"
set xr [-pi:pi]
set yr [-pi:pi]
set palette rgb 33,13,10

set size 2
set size square
unset key
set cblabel "kcal/mol"
set xlabel "{/Symbol f} (rad)"
set ylabel "{/Symbol y} (rad)"

p "femap.dat" w image, \
  "forces_gnuplot.dat" u 1:2 w p pt 3 lc rgbcolor "white", \
  "" u 1:2:($3/50):($4/50) w vectors lc rgbcolor "white", \
  "ref_md.contours" w l lc rgbcolor"black"
  

