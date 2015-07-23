#!/bin/bash
make
echo -n "PARAMETRI: < tipo_ds input_points dimensioni centroidi generazioni> "
read text
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
for i in 1 2 3
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
echo "#DEFC v9b"
time ./defc9b.x $text
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
done
