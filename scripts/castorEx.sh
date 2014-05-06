#!/bin/bash


Thr=MinBias
B=0

N=0

for i in `rfdir /castor/cern.ch/user/m/mmarionn/Ntuples/v3/ | awk '{print $9}'`; do

    if [ $i == $Thr ] && [ $B -eq 0 ]; then
	B=1
    fi


    if [ ! $B -eq 1 ]; then 
	continue
    fi


for j in `rfdir /castor/cern.ch/user/m/mmarionn/Ntuples/v3/$i |  awk '{print $9}'`; do


stager_get -M /castor/cern.ch/user/m/mmarionn/Ntuples/v3/$i/$j


N=`echo $N + 1 | bc`

if [ $N -eq 2950 ]; then 
    sleep 10m
N=0
fi


done

done
