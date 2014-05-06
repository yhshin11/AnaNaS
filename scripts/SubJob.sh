#!/bin/bash

BASE=/castor/cern.ch/user/m/mmarionn/
VERSION=ReprocessMET_v5 #/MC_v3 #_NoPUSmear_PUXSect 
FILTER=DoubleMu_178001_179000 #DoubleMu_163001_163869  #DoubleMu_178001_179000  #DoubleMu_163001_163869
N=0
OPATH=v18_2011_PhiMod  #mc_v14_2011


nohup ./CleanLSF.sh &

rfmkdir /castor/cern.ch/user/m/mmarionn/Ntuples/$OPATH

for i in `rfdir $BASE/$VERSION | awk '{print $9}'`; do


  #  REP=` rfdir /castor/cern.ch/user/m/mmarionn/Ntuples/$OPATH/$DATASET`

 


    DATASET=$i

 #   if [ "$DATASET" != "$FILTER" ]; then
#	continue
#    fi

   rfmkdir /castor/cern.ch/user/m/mmarionn/Ntuples/$OPATH/$i


    for j in `rfdir $BASE/$VERSION/$DATASET/ | awk '{print $9}'`; do

    FILE=$j
    MONTH=`rfdir $BASE/$VERSION/$DATASET/$j | awk '{print $6}'`
    DAY=`rfdir $BASE/$VERSION/$DATASET/$j | awk '{print $7}'`

    if [ "$N" == "300" ]; then
	break
    fi

    OFILE=` rfdir /castor/cern.ch/user/m/mmarionn/Ntuples/$OPATH/$DATASET/AnaTuple_$FILE`
#    echo " test    "$OFILE
    if [ "${OFILE:0:3}" == "-rw" ]; then
#	echo "boulou"
	continue
    fi
 #   echo "bili"

    bsub -R "pool>10000" -q 8nh source castor.sh $VERSION $DATASET $FILE Diboson $OPATH

    N=`echo $N + 1 | bc`

   
    #break
   

    done


done

