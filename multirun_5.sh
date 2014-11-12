#!/bin/bash
make
echo -n "PARAMETRI: < tipo_ds input_points dimensioni centroidi generazioni tipo_dw > "
read text
echo "PARAMETRI LETTI: $text"
for i in 1 2 3 4 5
do
echo "run # $i"
#echo "DEFC v5"
#time ./defc5 $text
#echo "DEFC v3"
#time ./defc3 $text
#echo "DEFC v2"
#time ./defc2 $text
echo "DEFC v7"
time ./defc7 $text
done
./plot.sh
