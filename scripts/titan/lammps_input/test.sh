#! /bin/bash

i=16
while [ $i -le 256 ]
do
    cp in.lj.32v16.txt in.lj.32v16_$i.txt
    i=$(( $i * 2 ))
    python test.py
done 
