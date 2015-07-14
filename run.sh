#!/bin/bash
echo "making..."
make
echo -n "PARAMETRI: < tipo_ds input_points dimensioni centroidi generazioni > "
read text
echo "#DEFC v9b"
time ./defc9b.x $text &
