#!/bin/bash
make
echo -n "PARAMETRI: < tipo_ds input_points dimensioni centroidi generazioni tipo_dw > "
read text
echo "letti: $text"
for i in 1 2
do
   echo "run # $i"
#./defc5 3 1000 2 5 1000 2
./defc5 $text
done
