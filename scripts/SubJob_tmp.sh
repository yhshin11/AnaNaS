#!/bin/bash

BASE=/castor/cern.ch/user/m/mmarionn/
VERSION=ReprocessMET_v2/MC_v4
FILTER=Z_2m_HCAL25
N=0

OUTPATH=

for i in `rfdir $BASE/$VERSION | awk '{print $9}'`; do


    rfmkdir /castor/cern.ch/user/m/mmarionn/Ntuples/v3/$i


    DATASET=$i

 #   if [ "$DATASET" != "$FILTER" ]; then
#	continue
#    fi


    for j in `rfdir $BASE/$VERSION/$DATASET/ | awk '{print $9}'`; do

    FILE=$j


    OFILE=` rfdir /castor/cern.ch/user/m/mmarionn/Ntuples/v3/$DATASET/AnaTuple_$FILE`
#    echo " test    "$OFILE
    if [ "${OFILE:0:3}" == "-rw" ]; then
#	echo "boulou"
	continue
    fi
 #   echo "bili"

    bsub -R "pool>10000" -q 1nh source castor_tmp.sh $VERSION $DATASET $FILE Diboson

    N=`echo $N + 1 | bc`

   
    #break
   

    done


done

