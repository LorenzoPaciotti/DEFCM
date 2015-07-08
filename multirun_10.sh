#!/bin/bash
make
echo -n "PARAMETRI: < tipo_ds input_points dimensioni centroidi generazioni tipo_dw > "
read text
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
for i in 1 2 3 4
do
echo "run # $i"
#echo "DEFC v2"
#time ./defc2 $text
#echo "DEFC v3"
#time ./defc3 $text
#echo "#DEFC v5"
#./defc5 $text
#echo "#DEFC v7"
#time ./defc7 $text
echo "#DEFC v9"
time ./defc9 $text
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
done
