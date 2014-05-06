#!/bin/bash

BASE=/eos/cms/store/user/mmarionn/AnaNaSNtuples/   #/Castor/cern.ch/user/m/mmarionn/
VERSION=METv2  #Prod8TeV_v4   #Prod8TeV_53X_v1
ListDS=DatasetFile

#ls $ListDS

ANTYPE=MET #MET

TMPOPATH=/afs/cern.ch/user/m/mmarionn/workspace/private/ANSNtuples53Xv4 #ANSNtuples53XABCDMu/   #ANSNtuples53Xv1/
OVERSION=.

OPATH=$TMPOPATH/$OVERSION

ANPATH=$ANANAS
N=0
nohup ./CleanLSF.sh &

if [ ! -d $OPATH ]; then
    mkdir -p $OPATH
fi


while read line; do


echo $line
  #  REP=` rfdir /castor/cern.ch/user/m/mmarionn/Ntuples/$OPATH/$DATASET`
    
    DATASET=$line
    
#FIXME
 #   if [ "$DATASET" == "TTbarJets" ]; then
#	VERSION=Prod8TeV_v2
#    fi


    if [ "${DATASET:0:2}" == "//" ]; then
	echo "skipping $DATASET"
	continue
    fi


    if [ -z "$OPATH/$line" ]; then
	mkdir $OPATH/$line
    fi

    echo $BASE/$VERSION/$DATASET/

#      for j in `rfdir $BASE/$VERSION/$DATASET/ | awk '{print $9}'`; do
    for j in `eos ls $BASE/$VERSION/$DATASET/`; do  #/afs/cern.ch/project/eos/installation/0.2.5/bin/eos.select

    	FILE=$j
#	MONTH=`rfdir $BASE/$VERSION/$DATASET/$j | awk '{print $6}'`
#	DAY=`rfdir $BASE/$VERSION/$DATASET/$j | awk '{print $7}'`
	
	OFILE=$OPATH/$DATASET/AnaTuple_$FILE

	#echo $OFILE

	#protection against existing files
	if [ -e $OFILE ]; then
	  echo "file existing : $OFILE" 
	  continue
	fi

#	echo $FILE
	
	#echo $VERSION $DATASET $FILE $ANTYPE $OPATH $ANPATH
	bsub -R "pool>10000" -q 1nh source singleJob.sh $BASE/$VERSION $DATASET $FILE $ANTYPE $OPATH $ANPATH
	# source singleJob.sh $BASE/$VERSION $DATASET $FILE $ANTYPE $OPATH $ANPATH

	N=`echo $N + 1 | bc`

	if [[ $N == 3000 ]]; then
	    N=0
	    sleep 10m
	fi
	
	#break;

    done


done < $ListDS

