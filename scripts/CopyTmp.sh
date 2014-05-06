#!/bin/bash

BASE=/castor/cern.ch/user/m/mmarionn/
VERSION=ReprocessMET_v2/MC_v4/
FILTER=Z_2m
N=0

OUTPATH=

for i in `rfdir $BASE/$VERSION | awk '{print $9}'`; do


  #  rfmkdir /castor/cern.ch/user/m/mmarionn/Ntuples/v5/$i
    mkdir /tmp/mmarionn/$i

    DATASET=$i

    if [ "${DATASET:0:6}" != "$FILTER" ]; then
	continue
    fi


    for j in `rfdir $BASE/$VERSION/$DATASET/ | awk '{print $9}'`; do

    FILE=$j
    MONTH=`rfdir $BASE/$VERSION/$DATASET/$j | awk '{print $6}'`
    DAY=`rfdir $BASE/$VERSION/$DATASET/$j | awk '{print $7}'`

    if [ "$MONTH" != "Dec" ]; then
	continue
    fi

    rfcp $BASE/$VERSION/$DATASET/$j /tmp/mmarionn/$DATASET/$j

 #   OFILE=` rfdir /castor/cern.ch/user/m/mmarionn/Ntuples/v5/$DATASET/AnaTuple_$FILE`
#    echo " test    "$OFILE
  #  if [ "${OFILE:0:3}" == "-rw" ]; then
#	echo "boulou"
#	continue
#    fi
 #   echo "bili"

#    bsub -R "pool>10000" -q 1nh source castor.sh $VERSION $DATASET $FILE Diboson

    N=`echo $N + 1 | bc`

   
    #break
   

    done


done

