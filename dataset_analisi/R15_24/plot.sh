#!/usr/bin/gnuplot -p
#plot 'x.dat'
set terminal postscript eps enhanced color "Helvetica" 20
set output "R15.eps"

 
# beginning of multiplot
set multiplot


set grid
set samples 10000
 
# the first plot
#set size 1,1
#set origin 0,0
#set xlabel "X"
#set ylabel "Y"
set xrange[0:20]
set yrange[0:20]
unset key
plot "x.dat" with p pt 1
unset grid

set key center top
plot "v_defc.out" with p pt 11 ps 3 lc 0

set key left top
plot "v_fcm.out" with p pt 15 ps 2 lc 15
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
