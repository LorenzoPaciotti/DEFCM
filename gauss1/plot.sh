#!/usr/bin/gnuplot
 
set terminal postscript eps enhanced color "Helvetica" 20
set output "gauss4.eps"

 
# beginning of multiplot
set multiplot


set grid
 
# the first plot
#set size 1,1
#set origin 0,0
#set xlabel "X"
#set ylabel "Y"
set xrange[-5:70]
set yrange[-5:70]
unset key
plot "x.dat" with p pt 1
unset grid

set key top center
plot "v_defc7.dat" with p pt 11 ps 3 lc 0

set key top left
plot "v_fcm.dat" with p pt 15 ps 2 lc 15
# smaller overlay plot
#set size 0.40,0.40
#set origin 0.5,0.5
#set xrange[-5:70]
#set yrange[-5:70]
#set lmargin 5
#unset xlabel
#unset ylabel
#unset arrow
#unset logscale
#set xtics 0.05
#set ytics 0.2
#set grid
#plot "v_defc.dat" title "matrice V" with p pt 3
 
# end of multiplot
#unset multiplot
