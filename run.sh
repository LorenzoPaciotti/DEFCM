#!/bin/bash
git pull
make
echo -n "PARAMETRI: < tipo_ds input_points dimensioni centroidi generazioni tipo_dw > "
read text
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
echo "#DEFC v9"
time ./defc9 $text &
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
echo "#DEFC v9b"
time ./defc9b $text &
