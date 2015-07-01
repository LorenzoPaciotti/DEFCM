#!/bin/bash
make
echo -n "PARAMETRI: < tipo_ds input_points dimensioni centroidi generazioni tipo_dw > "
read text
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
for i in 1 2 3 4 5 6 7 8 9 10
do
echo "run # $i"
#echo "#DEFC v5"
#./defc5 $text
#echo "DEFC v3"
#time ./defc3 $text
#echo "DEFC v2"
#time ./defc2 $text
#echo "#DEFC v7"
#time ./defc7 $text
time ./defc9 $text
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
done
