#!/usr/bin/gnuplot
 
set terminal postscript eps enhanced color "Helvetica" 20
set output 'my_data.eps'
 
# beginning of multiplot
set multiplot

set key off
 
# the first plot
set size 1,1
set origin 0,0
set xlabel "X"
set ylabel "Y"
set xrange[-5:70]
set yrange[-5:70]
plot "x.dat" title "matrice X" with p pt 5
 
# smaller overlay plot
set size 0.45,0.45
set origin 0.15,0.5
set xrange[-5:70]
set yrange[-5:70]
set lmargin 5
unset xlabel
unset ylabel
unset arrow
unset logscale
#set xtics 0.05
#set ytics 0.2
set grid
plot "v_defc.dat" title "matrice V" with p pt 5
 
# end of multiplot
unset multiplot
