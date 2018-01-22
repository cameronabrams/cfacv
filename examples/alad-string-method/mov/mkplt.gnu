#!/usr/bin/gnuplot
set term post eps enh color "Helvetica" 24 linewidth 2
set out "tmp.eps"
set size 1.5,1.5 
set size square

set xr [-pi:pi]
set yr [-pi:pi]
set zr [0:16]
set cbr [0:16]

set xlabel "{/Symbol f} (rad)"
set ylabel "{/Symbol y} (rad)" offset 0.5,0

unset key

set sty l 1 lc rgb "blue" pt 7 ps 1.5 lw 2
set sty l 2 lc rgb "yellow" pt 7 ps 1.85 lw 1

set format cb "%.1f"

set palette rgbformulae 33,13,10
set palette maxcolors 48 

s=0.1

p "../../../alad-single-sweep/output/femap.dat" not w image,\
  "s.dat" not w p ls 2, "s.dat" not w lp ls 1
#,\
#  "g.dat" u 1:2:(-$3*s):(-$4*s) not w vec lc rgbcolor "white"

#  "../../../alad-single-sweep/ref_md/ref_md.contours" not w l lc rgbcolor "black",\


