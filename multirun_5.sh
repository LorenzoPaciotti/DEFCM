#!/bin/bash
make
echo -n "PARAMETRI: < tipo_ds input_points dimensioni centroidi generazioni tipo_dw > "
read text
echo "letti: $text"
#for i in 1 2
#do
#   echo "run # $i"
#./defc5 3 1000 2 5 1000 2
echo "DEFC v5"
./defc5 $text
echo "DEFC v3"
./defc3 $text
echo "DEFC v2"
./defc2 $text
#done
./plot.sh
