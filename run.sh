#!/bin/bash
echo "controllo git"
git status
git pull
echo "making..."
make
echo -n "PARAMETRI: < tipo_ds input_points dimensioni centroidi generazioni > "
read text
#echo "#DEFC v9"
#time ./defc9 $text &
echo "#DEFC v9b"
time ./defc9b $text &
top
