#!/usr/bin/gnuplot

set terminal postscript enhanced color
set output '| ps2pdf - gauss5.pdf'
set multiplot
set grid
set xrange[-5:75]
set yrange[-5:75]
unset key
plot "x.dat" with p pt 1
unset grid

set key top center
plot "v_defc.dat" with p pt 11 ps 3 lc 0

set key top left
plot "v_fcm.dat" with p pt 15 ps 2 lc 15
