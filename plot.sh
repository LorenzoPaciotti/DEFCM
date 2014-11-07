#!/usr/bin/gnuplot
 
set terminal postscript eps enhanced color "Helvetica" 20

 
# beginning of multiplot
#set multiplot

set key off
set grid
 
# the first plot
set size 1,1
set origin 0,0
#set xlabel "X"
#set ylabel "Y"
set xrange[-10:50]
set yrange[-10:50]
set output 'X.eps'
plot "x.dat" title "matrice X" with p pt 3

set terminal postscript eps enhanced color "Helvetica" 20
set output 'V_DEFCM.eps'
plot "v_defc.dat" title "matrice X" with p pt 3 ps 3

set terminal postscript eps enhanced color "Helvetica" 20
set output 'V_FCM.eps'
plot "v_fcm.dat" title "matrice X" with p pt 3 ps 3
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
